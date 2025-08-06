export mpmsolver!
export procedure!

include(joinpath(@__DIR__, "solver/interpolation.jl"))
include(joinpath(@__DIR__, "solver/datatransfer.jl"))
include(joinpath(@__DIR__, "solver/kernels.jl"))
include(joinpath(@__DIR__, "solver/materials.jl"))
include(joinpath(@__DIR__, "solver/calculator.jl"))


mpmsolver!(solver::Function, args...) = solver(args...)

function procedure!(conf::Config, grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
    t_cur = T2(conf.t_cur)
    t_tol = T2(conf.t_tol)
    Δt    = T2(conf.Δt)
    dev   = conf.dev
    h5    = conf.h5
    dev_grid, dev_mpts = host2device(dev, grid, mpts)
    memsize = memorysize(dev_grid) + memorysize(dev_mpts)
    @info "uploaded $(@sprintf("%.2f", memsize)) GiB data to $(dev)"

    fid = set_hdf5(conf)
    printer = set_pb(conf)
    conf.stime[] = time()
    while t_cur < t_tol
        hdf5!(h5, fid, t_cur, mpts, dev_mpts)
        resetgridstatus!(dev_grid)
        resetmpstatus!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts)
        p2g!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts)
        solvegrid!(dev)(ndrange=dev_grid.ni, dev_grid, Δt)
        doublemapping1!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, Δt)
        doublemapping2!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts)
        doublemapping3!(dev)(ndrange=dev_grid.ni, dev_grid, Δt)
        g2p!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts)
        pb!(printer, t_cur, Δt)
        t_cur += Δt
        h5.iters[] += 1
    end
    conf.etime[] = time(); KAsync(dev)
    device2host!(mpts, dev_mpts)
    @info "downloaded $(@sprintf("%.2f", memorysize(dev_mpts))) GiB data to CPU()"
    hdf5!(h5, fid, grid)
    close(fid)
end