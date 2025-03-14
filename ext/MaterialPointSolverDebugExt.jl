module MaterialPointSolverDebugExt

using MaterialPointSolver
using KernelAbstractions
using WGLMakie

import MaterialPointSolver: materialpointsolver!, info_print, submit_work!, perf, procedure!
import MaterialPointSolver: DeviceArgs, DeviceGrid, DeviceParticle, DeviceProperty, 
                            DeviceVBoundary, DebugConfig
import MaterialPointSolver: initmpstatus!, progressinfo, host2device, updatepb!, KAsync,
                            getBackend, clean_device!
import StatsBase: mean

# custom dark theme for WGLMakie LScene
custom_dark() = Theme(
    Axis = (backgroundcolor="transparent",
            bottomspinevisible=false,
            leftspinevisible=false,
            rightspinevisible=false,
            topspinevisible=false,
            xgridcolor=(:white, 0.09),
            xlabelpadding=3,
            xminorticksvisible=false,
            xticksvisible=false,
            ygridcolor=(:white, 0.09),
            ylabelpadding=3,
            yminorticksvisible=false,
            yticksvisible=false),
    Axis3 = (xgridcolor=(:white),
             xspinesvisible=false,
             xticksvisible=false,
             ygridcolor=(:white),
             yspinesvisible=false,
             yticksvisible=false,
             zgridcolor=(:white),
             zspinesvisible=false,
             zticksvisible=false),
    backgroundcolor="#181818",
    Colorbar = (spinewidth=0,
                ticklabelpad=5,
                ticksvisible=false),
    Legend = (framevisible=false,
              padding=(0, 0, 0, 0)),
    linecolor=:white,
    textcolor=:white
)

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
    pdims = size(mp.ξ, 2)
    vtime = Observable(Ti)
    vpts  = Observable(Array{Float32, 2}(mp.ξ))
    vattr = Observable(Array{Float32, 1}(getproperty(mp, debug.plot.attr)))

    vmin = minimum(grid.ξ, dims=1)
    vmax = maximum(grid.ξ, dims=1)
    gridvertices = [
        Point3f(vmin[1], vmin[2], vmin[3]),
        Point3f(vmax[1], vmin[2], vmin[3]),
        Point3f(vmax[1], vmax[2], vmin[3]),
        Point3f(vmin[1], vmax[2], vmin[3]),
        Point3f(vmin[1], vmin[2], vmax[3]),
        Point3f(vmax[1], vmin[2], vmax[3]),
        Point3f(vmax[1], vmax[2], vmax[3]),
        Point3f(vmin[1], vmax[2], vmax[3])
    ]
    edges = [
        (1, 2), (2, 3), (3, 4), (4, 1),  # 底面
        (5, 6), (6, 7), (7, 8), (8, 5),  # 顶面
        (1, 5), (2, 6), (3, 7), (4, 8)   # 侧边
    ]
    line_vertices = Point3f[]
    for (i, j) in edges
        push!(line_vertices, gridvertices[i])
        push!(line_vertices, gridvertices[j])
    end

    set_theme!(custom_dark())
    fig = Figure()
    ax = LScene(fig[1, 1], show_axis=debug.plot.axis)
    Label(fig[1, 1, Top()], @lift("Simulation Time: $(lpad(round($vtime, digits=2), 8)) s"), padding = (0, 0, 5, 0))
    canvas = ax.scene.plots[1]
    canvas.ticks[:textcolor] = :white
    canvas.ticks[:fontsize] = (debug.plot.cbfontsize, 
                               debug.plot.cbfontsize, 
                               debug.plot.cbfontsize)
    canvas.frame[:axiscolor] = ("#818181", "#818181", "#818181")
    canvas.names[:textcolor] = (:white, :white, :white)
    canvas.names[:fontsize] = (debug.plot.axfontsize, 
                               debug.plot.axfontsize, 
                               debug.plot.axfontsize)
    linesegments!(ax, line_vertices, linewidth=0)
    plt = scatter!(ax, vpts, color=vattr, markersize=debug.plot.pointsize, 
        colormap=debug.plot.colormap, transparency=false, depthsorting=false, marker=:rect)
    if debug.plot.colorbar
        tickspace = debug.plot.tickspace
        tickprecision = debug.plot.tickprecision
        tickformat = debug.plot.tickformat
        Colorbar(fig[1, 2], plt, label=debug.plot.colorby, spinewidth=0, 
            tickformat="{:$(tickspace).$(tickprecision)$(tickformat)}")
    end

    center = mean(grid.ξ, dims=1)
    cam = cameracontrols(ax.scene)
    cam.lookat[] = pdims==2 ? Vec2f(center...) : Vec3f(center...)
    camoffset = debug.plot.camoffset
    cam.eyeposition[] = pdims==2 ? Vec2f(
        center[1] + camoffset, center[2] + camoffset) : Vec3f(
        center[1] + camoffset, center[2] + camoffset, center[3] + camoffset)
    cam.upvector[] = pdims==2 ? Vec2f(0, 1) : Vec3f(0, 0, 1)
    update_cam!(ax.scene, cam)
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
            copyto!(vpts[], dev_mp.ξ)
            copyto!(vattr[], getproperty(dev_mp, debug.plot.attr))
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