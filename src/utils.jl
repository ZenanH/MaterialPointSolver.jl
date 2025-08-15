
export format_seconds, set_pb, update_pb!, finish_pb!, set_hdf5, hdf5!, model_info

@inline function format_seconds(s_time)
    s = s_time < 1 ? 1.0 : ceil(Int, s_time)
    dt = Dates.Second(s)
    days = Dates.value(dt) ÷ (60 * 60 * 24)
    time = Dates.Time(Dates.unix2datetime(Dates.value(dt) % (60 * 60 * 24)))
    return days == 0 ?
        Dates.format(time, "HH:MM:SS") :
        @sprintf("%02d days: %s", days, Dates.format(time, "HH:MM:SS"))
end

@inline function set_pb(conf::Config)
    conf.stime[] = time()
    p = Progress(100; dt=conf.log_int,
        desc      = "\e[1;36m[ Info:\e[0m solving",
        barlen    = 20,
        barglyphs = BarGlyphs(" ■■  "),
        output    = stderr,
        enabled   = true
    )
    return p 
end

@inline function update_pb!(printer::Progress, t_cur, t_tol)
    pos = unsafe_trunc(Int, (t_cur / t_tol) * 100)
    pos < 100 && update!(printer, pos)
end

@inline function finish_pb!(conf::Config, printer::Progress)
    conf.etime[] = time()
    finish!(printer)
end

set_hdf5(conf::Config) = HDF5.h5open(joinpath(conf.prjdst, "$(conf.prjname).h5"), "w")

@inline _get_nested_field(obj, path) = foldl(getfield, path; init=obj)

@inline function hdf5!(conf::H5_T, fid, t_cur, mpts, dev_mpts)
    if conf.iters[] == 0 || conf.iters[] % conf.interval == 0
        device2host!(mpts, dev_mpts, conf.varnames)
        g = create_group(fid, "group$(conf.gname[])")
        @inbounds for path in conf.fpvar
            g[string(path[end])] = _get_nested_field(mpts, path)
        end
        g["time"] = t_cur
        conf.gname[] += 1
    end
end

@inline function hdf5!(::H5_T, fid, grid)
    g = create_group(fid, "grid")
    @inbounds for vars in (:h, :x1, :x2, :y1, :y2, :z1, :z2)
        g[string(vars)] = getfield(grid, vars)
    end
end

@inline hdf5!(::H5_F, fid, t_cur, mpts::DeviceParticle{T1, T2}, dev_mpts::DeviceParticle{T1, T2}) where {T1, T2} = nothing
@inline hdf5!(::H5_F, fid, grid::DeviceGrid{T1, T2}) where {T1, T2} = nothing


function model_info(conf::Config, grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
    textplace1 = 16
    textplace2 = 9
    basis = lpad(string(conf.basis), textplace1)
    h5 = typeof(conf.h5) <: H5_T ? lpad("true", textplace1) : lpad("false", textplace1)
    material = lpad(string(conf.material), textplace1)
    t_tol = rpad(string(@sprintf("%.2e", conf.t_tol))*" s", textplace2)
    npts = rpad(string(@sprintf("%.2e", mpts.np)), textplace2)
    ngid = rpad(string(@sprintf("%.2e", grid.ni)), textplace2)
    
    FLIP = lpad(string(@sprintf("%.2f", mpts.FLIP)), 6)
    PIC  = lpad(string(@sprintf("%.2f", 1-mpts.FLIP)), 6)
    precision = typeof(T2) == Float32 ? "single" : "double"

    
    @info """$(conf.prjname) [$(conf.dev)]
    ─────────────┬────────────────────────────┬──────────────────
    FLIP: $(FLIP) │ basis   : $(basis) │ mpts : $(npts)
    PIC : $(PIC) │ HDF5    : $(h5) │ nodes: $(ngid)
    ϵ   : $(precision) │ material: $(material) │ t_tol: $(t_tol)
    ─────────────┴────────────────────────────┴──────────────────
    """
end