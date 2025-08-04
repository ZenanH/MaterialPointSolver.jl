function launchkernels!(conf, dev_grid, dev_mpts, Δt)
    dev   = conf.dev
    basis = conf.basis
    dim   = conf.dim
    ϵ     = conf.ϵ
    G     = dev_mpts.G
    FLIP  = dev_mpts.FLIP
    resetgridstatus!(dev_grid, ϵ)
    resetmpstatus!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, basis, dim, ϵ)
    P2G!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, G, dim, ϵ)
    solvegrid!(dev)(ndrange=dev_grid.ni, dev_grid, Δt, dim, ϵ)
    doublemapping1!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, Δt, FLIP, dim, ϵ)
    doublemapping2!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, dim, ϵ)
    doublemapping3!(dev)(ndrange=dev_grid.ni, dev_grid, dim, ϵ)
    G2P!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, Δt, dim, ϵ)
end

function procedure!(conf, grid, mpts)
    t_cur = conf.t_cur
    t_tol = conf.t_tol
    Δt = conf.Δt
    dev = conf.dev
    h5 = conf.h5
    dev_grid, dev_mpts = host2device(dev, grid, mpts)
    memsize = totalsize(dev_grid, dev_mpts)
    @info "uploaded $(@sprintf("%.2f", memsize)) GiB data to $(dev)"

    fid = set_hdf5(conf)
    printer = set_pb(conf)
    conf.time_start[] = time()
    while t_cur < t_tol
        # hdf5 data export
        hdf5!(h5, fid, t_cur, mpts, dev_mpts)
        # mpm kernel functions
        launchkernels!(conf, dev_grid, dev_mpts, Δt)
        # terminal printer
        pb!(printer, t_cur, Δt)
        # update vars
        t_cur += Δt
        h5.iters[] += 1
    end
    KAsync(dev); conf.time_end[] = time()
    device2host!(mpts, dev_mpts)
    @info "downloaded $(@sprintf("%.2f", totalsize(dev_mpts))) GiB data to CPU()"
    hdf5!(conf.h5, fid, grid)
    close(fid)
end