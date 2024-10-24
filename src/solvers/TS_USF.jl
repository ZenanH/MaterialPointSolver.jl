#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : TS_USF.jl                                                                  |
|  Description: Defaule USF (update stress first) implementation for two-phase             |
|               single-point MPM                                                           |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Functions  : procedure! [2D & 3D]                                                       |
+==========================================================================================#

function procedure!(
    args::     DeviceArgs{T1, T2}, 
    grid::     DeviceGrid{T1, T2}, 
    mp  :: DeviceParticle{T1, T2},
    attr:: DeviceProperty{T1, T2}, 
    bc  ::DeviceVBoundary{T1, T2},
    ΔT  ::T2,
    Ti  ::T2,
        ::Val{:TS},
        ::Val{:USF}
) where {T1, T2}
    G = Ti < args.Te ? args.gravity / args.Te * Ti : args.gravity
    dev = getBackend(Val(args.device))
    G2P_TS!(dev)(ndrange=mp.np, grid, mp, ΔT)
    if args.constitutive == :hyperelastic
        hyE!(dev)(ndrange=mp.np, mp, attr)
    elseif args.constitutive == :linearelastic
        liE!(dev)(ndrange=mp.np, mp, attr)
    elseif args.constitutive == :druckerprager
        liE!(dev)(ndrange=mp.np, mp, attr)
        if Ti ≥ args.Te
            dpP!(dev)(ndrange=mp.np, mp, attr)
        end
    elseif args.constitutive == :mohrcoulomb
        liE!(dev)(ndrange=mp.np, mp, attr)
        if Ti ≥ args.Te
            mcP!(dev)(ndrange=mp.np, mp, attr)
        end
    end
    resetgridstatus_TS!(dev)(ndrange=grid.ni, grid)
    resetmpstatus_TS!(dev)(ndrange=mp.np, grid, mp, Val(args.basis))
    P2G_TS!(dev)(ndrange=mp.np, grid, mp, G)
    solvegrid_USL_TS!(dev)(ndrange=grid.ni, grid, bc, ΔT, args.ζs)
    doublemapping1_TS!(dev)(ndrange=mp.np, grid, mp, ΔT, args.FLIP, args.PIC)
    if args.MVL == true
        vollock1_TS!(dev)(ndrange=mp.np, grid, mp)
        vollock2_TS!(dev)(ndrange=mp.np, grid, mp)
    end                                  
    return nothing
end