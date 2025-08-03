#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: utils file                                                                 |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, ISTE, Université de Lausanne                                   |
|  Maintainer : Zenan Huo                                                                  |
+==========================================================================================#

export format_seconds, set_pb, pb!, set_hdf5, hdf5!

versioninfo() = println(
    "v" * string(Pkg.dependencies()[Base.UUID("80d0c2dc-c1ec-4aa7-9846-c2b398d9a94e")].version))

@inline function format_seconds(s_time)
    s = s_time < 1 ? 1.0 : ceil(Int, s_time)
    dt = Dates.Second(s)
    days = Dates.value(dt) ÷ (60 * 60 * 24)
    time = Dates.Time(Dates.unix2datetime(Dates.value(dt) % (60 * 60 * 24)))
    return days == 0 ?
        Dates.format(time, "HH:MM:SS") :
        @sprintf("%02d days: %s", days, Dates.format(time, "HH:MM:SS"))
end

set_pb(conf) = (time1=Ref(time()), t_tol=conf.t_tol, interval=conf.log_interval, iters=Ref{Int}(0))

@inline function pb!(printer::NamedTuple, t_cur, Δt)
    time2 = time()
    t_dur = time2 - printer.time1[]
    if t_dur > printer.interval
        speed   = max(printer.iters[] / t_dur, 1)
        eta_str = format_seconds(((printer.t_tol - t_cur) / Δt) / speed)
        time_str  = @sprintf("%.2f s", t_cur)
        speed_str = @sprintf("%.2f iters/s", speed)
        invo_str = "   \e[1;32m⇌\e[0m   "
        contents = "tcur: " * time_str * "$(invo_str)speed: " * speed_str * "$(invo_str)eta: " * eta_str
        print(stdout, "\r\e[2K"); print(stdout, "\e[1;32m[ Info: \e[0m", contents)
        printer.time1[] = time()
        printer.iters[] = 0
    end
    printer.iters[] += 1
    (t_cur + Δt ≥ printer.t_tol) && println()          # 只在结束时换行
end

set_hdf5(conf) = HDF5.h5open(joinpath(conf.prjdst, "$(conf.prjname).h5"), "w")

@inline function hdf5!(conf::hasHDF5, fid, t_cur, mpts, dev_mpts)
    if conf.iters[] == 0 || conf.iters[] % conf.interval == 0
        device2host!(mpts, dev_mpts, var=conf.varnames)
        g = create_group(fid, "group$(conf.gname[])")
        @inbounds for vars in conf.varnames
            g[string(vars)] = mpts[vars]
        end
        g["time"] = t_cur
        conf.gname[] += 1
    end
end

@inline function hdf5!(conf::hasHDF5, fid, grid)
    g = create_group(fid, "grid")
    fields = :z1 in keys(grid) ? (:h, :x1, :x2, :y1, :y2, :z1, :z2) : (:h, :x1, :x2, :y1, :y2)
    @inbounds for vars in fields
        g[string(vars)] = grid[vars]
    end
end

@inline hdf5!(conf::noHDF5, fid, t_cur, mpts, dev_mpts) = nothing
@inline hdf5!(conf::noHDF5, fid, grid) = nothing