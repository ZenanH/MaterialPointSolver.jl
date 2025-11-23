module Mechanics

using ..MaterialPointSolver

include(joinpath(@__DIR__, "funcs.jl"))

export mpm_mechanics!

function mpm_mechanics!(conf::Config, grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
    model_info(conf, grid, mpts)
    t_cur = T2(conf.t_cur)
    t_tol = T2(conf.t_tol)
    t_eld = T2(conf.t_eld)
    Δt    = T2(conf.Δt)
    dev   = conf.dev
    h5    = conf.h5
    dev_grid, dev_mpts = host2device(dev, grid, mpts)

    fid = set_hdf5(conf, mpts)
    printer = set_pb(conf)
    while t_cur < t_tol
        # main MPM loop
        Gg = 0 < t_eld < t_cur ? (mpts.G * t_cur) / t_eld : mpts.G
        resetgridstatus!(dev_grid)
        p2g!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, Gg)
        solvegrid!(dev)(ndrange=dev_grid.ni, dev_grid, Δt)
        doublemapping1!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts)
        doublemapping2!(dev)(ndrange=dev_grid.ni, dev_grid, Δt)
        g2p!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, t_eld, t_cur, Δt)

        # safty check for the simulation
        # status_checker(grid, mpts)
        
        # adaptive time step or fixed time step
        # Δt = conf.αT * reduce(min, dev_mpts.cfl) # adaptive time step
        Δt = conf.Δt # fixed time step

        # HDF5 output
        Δt = hdf5!(h5, fid, t_cur, t_tol, Δt, mpts, dev_mpts)
        
        t_cur += Δt
        update_pb!(printer, t_cur, t_tol)
    end
    finish_pb!(conf, printer); KAsync(dev)
    device2host!(mpts, dev_mpts)
    hdf5!(h5, fid, grid)
    close(fid)
end

end