export format_seconds
export set_pb, update_pb!, finish_pb!
export set_hdf5, hdf5, hdf5_finalizer
export model_info
export status_checker
export exit_sim
export start_sim
export timestep

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
        barglyphs = BarGlyphs(" ━━  "), # ■■
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

@inline _get_nested_field(obj, path) = foldl(getfield, path; init=obj)

function set_hdf5(conf::Config, mpts::DeviceParticle{T1, T2}) where {T1, T2}
    fid = HDF5.h5open(joinpath(conf.prjdst, "$(conf.prjname).h5"), "w")
    if typeof(conf.h5) <: H5_T
        g = create_group(fid, "group$(conf.h5.gname[])")
        @inbounds for path in conf.h5.fpvar
            g[string(path[end])] = _get_nested_field(mpts, path)
        end
        g["time"] = conf.t_cur
        conf.h5.gname[] += 1
        conf.h5.k[] += 1
    end
    return fid
end

@inline function hdf5(h5::H5_T, fid, t_cur, t_tol, Δti, mpts, dev_mpts)
    # check whether the current time has reached or surpassed the next save point.
    if h5.k[] ≤ h5.tol_iters && t_cur ≥ h5.interval[h5.k[]] - 1e-12
        device2host!(mpts, dev_mpts, h5.varnames)
        g = create_group(fid, "group$(h5.gname[])")
        @inbounds for path in h5.fpvar
            g[string(path[end])] = _get_nested_field(mpts, path)
        end
        g["time"] = t_cur
        h5.gname[] += 1
        h5.k[] += 1
    end
    
    # adjust the next time step to ensure it can precisely reach the next save point
    if h5.k[] ≤ h5.tol_iters
        next_save_time = h5.interval[h5.k[]]
        if t_cur < next_save_time
            # if the next step would exceed the save point, adjust Δt to exactly reach the save point
            if t_cur + Δti > next_save_time
                Δt = next_save_time - t_cur
            else
                Δt = Δti
            end
        end
    end
    
    # ensure not exceeding total time and Δt is valid
    Δt = min(Δt, t_tol - t_cur)
    
    return Δt
end

@inline function hdf5(h5::H5_T, fid, t_cur, t_tol, Δt_cfl, Δti, αT, mpts, dev_mpts)
    # correct wrong Δt_cfl
    if 0 < Δt_cfl < Δti
        Δt = clamp(Δt_cfl, Δti*αT, Δti) # using the cfl condition to determine the time step
    else
        Δt = αT*Δti    # using a fixed time step (including exceptional cases)
    end

    # check whether the current time has reached or surpassed the next save point.
    if h5.k[] ≤ h5.tol_iters && t_cur ≥ h5.interval[h5.k[]] - 1e-12
        device2host!(mpts, dev_mpts, h5.varnames)
        g = create_group(fid, "group$(h5.gname[])")
        @inbounds for path in h5.fpvar
            g[string(path[end])] = _get_nested_field(mpts, path)
        end
        g["time"] = t_cur
        h5.gname[] += 1
        h5.k[] += 1
    end
    
    # adjust the next time step to ensure it can precisely reach the next save point
    if h5.k[] ≤ h5.tol_iters
        next_save_time = h5.interval[h5.k[]]
        if t_cur < next_save_time
            # if the next step would exceed the save point, adjust Δt to exactly reach the save point
            if t_cur + Δt > next_save_time
                Δt = next_save_time - t_cur
            end
        end
    end
    
    # ensure not exceeding total time and Δt is valid
    Δt = min(Δt, t_tol - t_cur)
    
    return Δt
end

@inline function hdf5_finalizer(::H5_T, fid, grid)
    g = create_group(fid, "grid")
    @inbounds for vars in (:h, :x1, :x2, :y1, :y2, :z1, :z2)
        g[string(vars)] = getfield(grid, vars)
    end
end

@inline hdf5(::H5_F, fid, t_cur, mpts::DeviceParticle{T1, T2}, dev_mpts::DeviceParticle{T1, T2}) where {T1, T2} = nothing
@inline hdf5_finalizer(::H5_F, fid, grid::DeviceGrid{T1, T2}) where {T1, T2} = nothing
@inline hdf5(::H5_F, fid, t_cur, t_tol, Δti, mpts, dev_mpts) = Δti
@inline function hdf5(::H5_F, fid, t_cur, t_tol, Δt_cfl, Δti, αT, mpts, dev_mpts)
    # correct wrong Δt_cfl
    if 0 < Δt_cfl < Δti
        Δt = clamp(Δt_cfl, Δti*αT, Δti) # using the cfl condition to determine the time step
    else
        Δt = αT*Δti    # using a fixed time step (including exceptional cases)
    end
    
    return Δt
end


function model_info(conf::Config, grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
    textplace1 = 5
    FLIP = lpad(string(@sprintf("%.2f", mpts.FLIP)), textplace1)
    ζs   = lpad(string(@sprintf("%.2f", grid.ζs)), textplace1)
    h5   = typeof(conf.h5) <: H5_T ? lpad("true", textplace1) : lpad("false", textplace1)

    textplace2 = 9
    t_tol = rpad(string(@sprintf("%.2e", conf.t_tol))*" s", textplace2)
    t_eld = rpad(string(@sprintf("%.2e", conf.t_eld))*" s", textplace2)
    Δt = conf.adaptive ? rpad("adaptive  ", textplace2) :
        rpad(string(@sprintf("%.2e", conf.Δt))*" s", textplace2)

    textplace3 = 9
    npts = rpad(string(@sprintf("%.2e", mpts.np)), textplace3)
    ngid = rpad(string(@sprintf("%.2e", grid.ni)), textplace3)
    precision = typeof(T2) == Float32 ? rpad("single", textplace3) : rpad("double", textplace3)

    @info """\e[1;31mCommunity Edition\e[0m
    ────────────┬───────────────────┬────────────────
    FLIP: $(FLIP) │ t_tol: $(t_tol) │ mpts: $(npts)
      ζs: $(ζs) │ t_eld: $(t_eld) │ node: $(ngid)
    HDF5: $(h5) │    Δt: $(Δt) │    ϵ: $(precision)
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

function exit_sim(
    conf::Config, 
    printer, fid,
    grid::DeviceGrid{T1, T2}, 
    mpts::DeviceParticle{T1, T2}, 
    dev_mpts::DeviceParticle{T1, T2}
) where {T1, T2}
    finish_pb!(conf, printer)
    KAsync(conf.dev)
    device2host!(mpts, dev_mpts)
    hdf5_finalizer(conf.h5, fid, grid)
    close(fid)

    @info """simulation completed:
    $(conf.prjdst)/$(conf.prjname)
    """
end

function start_sim(
    conf::Config, 
    grid::DeviceGrid{T1, T2}, 
    mpts::DeviceParticle{T1, T2}
) where {T1, T2}
    model_info(conf, grid, mpts)
    t_cur = T2(conf.t_cur)
    t_tol = T2(conf.t_tol)
    t_eld = T2(conf.t_eld)
    Δt    = T2(conf.Δt)
    dev = conf.dev
    dev_grid, dev_mpts = host2device(dev, grid, mpts)

    fid = set_hdf5(conf, mpts)
    printer = set_pb(conf)

    return t_cur, t_tol, t_eld, Δt, dev, dev_grid, dev_mpts, fid, printer
end

function timestep(
    conf::Config, 
    mpts::DeviceParticle{T1, T2},
    dev_mpts::DeviceParticle{T1, T2},
    fid,
    t_cur::T2,
    t_tol::T2,
) where {T1, T2}
    # safty check for the simulation
    # status_checker(grid, mpts)
    if conf.adaptive
        Δt_cfl = conf.αT * reduce(min, dev_mpts.cfl)
        Δt = hdf5(conf.h5, fid, t_cur, t_tol, Δt_cfl, conf.Δt, conf.αT, mpts, dev_mpts)
    else
        Δt = hdf5(conf.h5, fid, t_cur, t_tol, conf.Δt, mpts, dev_mpts)
    end
    return Δt
end