module Mechanics

using ..MaterialPointSolver

include(joinpath(@__DIR__, "funcs.jl"))

export mpm_mechanics!

function mpm_mechanics!(conf::Config, grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
    t_cur, t_tol, t_eld, Δt, dev, dev_grid, dev_mpts, fid, printer = start_sim(conf, grid, mpts)
    while t_cur < t_tol # main MPM loop
        Gg = 0 < t_eld < t_cur ? (mpts.G * t_cur) / t_eld : mpts.G
        resetgridstatus!(dev_grid)
        p2g!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, Gg)
        solvegrid!(dev)(ndrange=dev_grid.ni, dev_grid, Δt)
        doublemapping1!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts)
        doublemapping2!(dev)(ndrange=dev_grid.ni, dev_grid, Δt)
        g2p!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, t_eld, t_cur, Δt)
        
        Δt = timestep(conf, mpts, dev_mpts, fid, t_cur, t_tol)
        t_cur += Δt
        update_pb!(printer, t_cur, t_tol)
    end
    exit_sim(conf, printer, fid, grid, mpts, dev_mpts)
end

end