#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : OS_MUSL.jl                                                                 |
|  Description: Defaule MUSL (modified update stress last) implementation for one-phase    |
|               single-point MPM                                                           |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
+==========================================================================================#

function procedure!(
    args::     DeviceArgs{T1, T2}, 
    grid::     DeviceGrid{T1, T2}, 
    mp  :: DeviceParticle{T1, T2}, 
    attr:: DeviceProperty{T1, T2},
    bc  ::DeviceVBoundary{T1, T2},
    ΔT  ::T2,
    Ti  ::T2,
        ::Val{:OS},
        ::Val{:MUSL}
) where {T1, T2}
    G = Ti < args.Te ? args.gravity / args.Te * Ti : args.gravity
    dev = getBackend(Val(args.device))
    # MPM procedure
    resetgridstatus_OS!(dev)(ndrange=grid.ni, grid)
    args.device == :CUDA && args.basis == :uGIMP ?
        resetmpstatus_OS_CUDA!(dev)(ndrange=mp.np, grid, mp, Val(args.basis)) :
        resetmpstatus_OS!(dev)(ndrange=mp.np, grid, mp, Val(args.basis))
    P2G_OS!(dev)(ndrange=mp.np, grid, mp, G)
    solvegrid_OS!(dev)(ndrange=grid.ni, grid, bc, ΔT, args.ζs)
    doublemapping1_OS!(dev)(ndrange=mp.np, grid, mp, attr, ΔT, args.FLIP, args.PIC)
    doublemapping2_OS!(dev)(ndrange=mp.np, grid, mp)
    doublemapping3_OS!(dev)(ndrange=grid.ni, grid, bc, ΔT)
    # F-bar based volumetric locking elimination approach
    if args.MVL == false
        G2P_OS!(dev)(ndrange=mp.np, grid, mp, ΔT)
    else
        G2Pvl1_OS!(dev)(ndrange=mp.np, grid, mp)
        fastdiv!(dev)(ndrange=grid.ni, grid)
        G2Pvl2_OS!(dev)(ndrange=mp.np, grid, mp, ΔT)
    end
    # update stress status
    if args.constitutive == :hyperelastic
        hyE!(dev)(ndrange=mp.np, mp, attr)
    elseif args.constitutive == :linearelastic
        liE!(dev)(ndrange=mp.np, mp, attr)
    elseif args.constitutive == :druckerprager
        liE!(dev)(ndrange=mp.np, mp, attr)
        Ti ≥ args.Te && dpP!(dev)(ndrange=mp.np, mp, attr)
    elseif args.constitutive == :mohrcoulomb
        liE!(dev)(ndrange=mp.np, mp, attr)
        Ti ≥ args.Te && mcP!(dev)(ndrange=mp.np, mp, attr)
    elseif args.constitutive == :bingham
        Ti < args.Te && liE!(dev)(ndrange=mp.np, mp, attr)
        Ti ≥ args.Te && bhP!(dev)(ndrange=mp.np, mp, attr, inv(ΔT))
    end
    return nothing
end