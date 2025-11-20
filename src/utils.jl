export format_seconds, set_pb, update_pb!, finish_pb!, set_hdf5, hdf5!, model_info, status_checker

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
        enabled   = true,
        color     = :cyan
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
    textplace1 = 5
    FLIP = lpad(string(@sprintf("%.2f", mpts.FLIP)), textplace1)
    ζs   = lpad(string(@sprintf("%.2f", grid.ζs)), textplace1)
    h5   = typeof(conf.h5) <: H5_T ? lpad("true", textplace1) : lpad("false", textplace3)

    textplace2 = 9
    t_tol = rpad(string(@sprintf("%.2e", conf.t_tol))*" s", textplace2)
    t_eld = rpad(string(@sprintf("%.2e", conf.t_eld))*" s", textplace2)
    Δt    = rpad(string(@sprintf("%.2e", conf.Δt))*" s", textplace2)

    textplace3 = 9
    npts = rpad(string(@sprintf("%.2e", mpts.np)), textplace3)
    ngid = rpad(string(@sprintf("%.2e", grid.ni)), textplace3)
    precision = typeof(T2) == Float32 ? rpad("single", textplace3) : rpad("double", textplace3)

    @info """\e[1;31mCommunity Edition\e[0m
    ────────────┬───────────────────┬────────────────
    FLIP: $(FLIP) │ t_tol: $(t_tol) │ mpts: $(npts)
    ζs  : $(ζs) │ t_eld: $(t_eld) │ node: $(ngid)
    HDF5: $(h5) │ Δt   : $(Δt) │ ϵ   : $(precision)
    ────────────┴───────────────────┴────────────────
    Project name: $(conf.prjname)
    """
end

@views @inbounds function status_checker(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}

    # 1. Ω check: cannot be negative
    if reduce(min, mpts.Ω) < 0
        error("MPM instability: Ω contains negative value")
    end

    # 2. ξ check: cannot contain NaN or Inf
    minξ = reduce(min, mpts.ξ)
    maxξ = reduce(max, mpts.ξ)
    if isnan(minξ) || isnan(maxξ) || isinf(minξ) || isinf(maxξ)
        error("MPM instability: ξ contains NaN or Inf")
    end

    # 3. σij check: cannot contain NaN or Inf
    minσ = reduce(min, mpts.σij)
    maxσ = reduce(max, mpts.σij)
    if isnan(minσ) || isnan(maxσ) || isinf(minσ) || isinf(maxσ)
        error("MPM instability: σij contains NaN or Inf")
    end

    # 4. particle position check: cannot exceed grid boundary
    mpts.vmin .= Array(reduce(min, mpts.ξ, dims=1))
    mpts.vmax .= Array(reduce(max, mpts.ξ, dims=1))
    
    rst = mpts.vmin[1] ≥ grid.x1 + grid.h * 2 &&
          mpts.vmax[1] ≤ grid.x2 - grid.h * 2 &&
          mpts.vmin[2] ≥ grid.y1 + grid.h * 2 &&
          mpts.vmax[2] ≤ grid.y2 - grid.h * 2 &&
          mpts.vmin[3] ≥ grid.z1 + grid.h * 2 &&
          mpts.vmax[3] ≤ grid.z2 - grid.h * 2
    if !rst error("MPM instability: particle position exceeds grid boundary") end

    return nothing
end