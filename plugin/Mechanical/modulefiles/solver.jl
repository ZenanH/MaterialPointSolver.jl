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