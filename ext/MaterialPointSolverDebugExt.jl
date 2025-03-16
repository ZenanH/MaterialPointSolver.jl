module MaterialPointSolverDebugExt

using MaterialPointSolver
using KernelAbstractions
using WGLMakie

import MaterialPointSolver: materialpointsolver!, info_print, submit_work!, perf, procedure!
import MaterialPointSolver: DeviceArgs, DeviceGrid, DeviceParticle, DeviceProperty, 
                            DeviceVBoundary, DebugConfig
import MaterialPointSolver: initmpstatus!, progressinfo, host2device, updatepb!, KAsync,
                            getBackend, clean_device!
import StatsBase: mean, sample

include(joinpath(@__DIR__, "DebugExt/plotutils.jl"))

"""
    materialpointsolver!(args::DeviceArgs, grid::DeviceGrid, mp::DeviceParticle, 
        attr::DeviceProperty, bc::DeviceVBoundary, debug::DebugConfig; 
        workflow::F=procedure!)

Description:
---
This function is the main entry point for the material point method (MPM) solver. Here we 
added `debug` to enable an in-situ visulaization of the simulation. This can help to check
the problem of the code. Note that it is just used for debug not for production mode. 
"""
function materialpointsolver!(
    args    ::     DeviceArgs{T1, T2}, 
    grid    ::     DeviceGrid{T1, T2}, 
    mp      :: DeviceParticle{T1, T2}, 
    attr    :: DeviceProperty{T1, T2},
    bc      ::DeviceVBoundary{T1, T2},
    debug   ::DebugConfig;
    workflow::F=procedure!
) where {T1, T2, F<:Function}
    info_print(args, grid, mp) # terminal info
    submit_work!(args, grid, mp, attr, bc, workflow, debug) # MPM solver
    perf(args) # performance summary
    return nothing
end

@views function submit_work!(
    args    ::     DeviceArgs{T1, T2},
    grid    ::     DeviceGrid{T1, T2}, 
    mp      :: DeviceParticle{T1, T2}, 
    attr    :: DeviceProperty{T1, T2},
    bc      ::DeviceVBoundary{T1, T2},
    workflow::F,
    debug   ::DebugConfig
) where {T1, T2, F<:Function}
    initmpstatus!(CPU())(ndrange=mp.np, grid, mp, Val(args.basis))
    # variables setup for the simulation 
    Ti = T2(0.0)
    pc = Ref{T1}(0)
    pb = progressinfo(args, "solving")
    ΔT = args.ΔT
    #ΔT = args.time_step==:auto ? cfl(args, grid, mp, attr, Val(args.coupling)) : args.ΔT
    dev_grid, dev_mp, dev_attr, dev_bc = host2device(grid, mp, attr, bc, Val(args.device))

    # plot setup
    vid = mp.np > debug.plot.pointnum ? sample(1:mp.np, debug.plot.pointnum, replace=false) : (1:mp.np)
    getvattr = debug.plot.calculate(grid, mp, attr, bc, vid)
    pdims = size(mp.ξ, 2)
    vtime = Observable(Ti)
    vpts  = Observable(Array{Float32, 2}(mp.ξ[vid, :]))
    vattr = Observable(Array{Float32, 1}(getvattr))

    set_theme!(gettheme())
    fig = Figure()
    ax = LScene(fig[1, 1], show_axis=debug.plot.axis)
    if debug.plot.axis
        canvas = ax.scene.plots[1]
        canvas.ticks[:textcolor] = :white
        canvas.ticks[:fontsize ] = debug.plot.axfontsize
        canvas.frame[:axiscolor] = "#818181"
        canvas.names[:textcolor] = :white
        canvas.names[:fontsize ] = debug.plot.axfontsize
    end
    Label(fig[1, 1, Top()], @lift("Simulation Time: $(lpad(round($vtime, digits=2), 8)) s"),
        padding = (0, 0, 5, 0))
    Label(fig[1, 1, Bottom()], "Debug mode: ON")
    linesegments!(ax, getbox(grid), linewidth=0)
    plt = scatter!(ax, vpts, color=vattr, markersize=debug.plot.pointsize, 
        colormap=debug.plot.colormap, transparency=false, depthsorting=false, marker=:rect)
    if debug.plot.colorbar
        Colorbar(fig[1, 2], plt, label=debug.plot.cbname, spinewidth=0, 
            tickformat=debug.plot.tickformat)
    end
    pdims==2 && cam2d!(ax.scene)
    display(fig)

    # main part: debug ON
    debug_id     = T1(1) # debug index
    debug_switch = T1(0) # debug step
    args.start_time = time()
    while Ti < args.Ttol
        workflow(args, dev_grid, dev_mp, dev_attr, dev_bc, ΔT, Ti, 
            Val(args.coupling), Val(args.scheme))
        ΔT = args.time_step == :auto ? args.αT * reduce(min, dev_mp.cfl) : args.ΔT
        Ti += ΔT
        debug_switch += 1
        args.iter_num += 1
        updatepb!(pc, Ti, args.Ttol, pb)
        if (debug_switch == args.hdf5_step) || (debug_switch == T1(0))
            copyto!(vpts[], dev_mp.ξ[vid, :])
            copyto!(vattr[], debug.plot.calculate(dev_grid, dev_mp, dev_attr, dev_bc, vid))
            notify(vpts)
            notify(vattr)
            vtime[] = Ti
            yield()
            debug_switch = T1(0)
            debug_id += T1(1)
        end
    end
    args.end_time = time()
    
    KAsync(getBackend(Val(args.device)))
    device2host!(grid, mp, attr, bc, dev_grid, dev_mp, dev_attr, dev_bc, Val(args.device))
    clean_device!(dev_grid, dev_mp, dev_attr, dev_bc, Val(args.device))
    return nothing
end

end