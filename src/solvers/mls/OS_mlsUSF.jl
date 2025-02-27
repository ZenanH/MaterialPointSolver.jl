#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : OS_mlsUSF.jl                                                               |
|  Description: Efficient MLS-MPM procedure                                                |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Notes      : Only quadratic basis function can be used in MLS-MPM, and we don't         |
|               calculate the gradient of the basis function                               |
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
        ::Val{:MLS}
) where {T1, T2}
    G = Ti < args.Te ? args.gravity / args.Te * Ti : args.gravity
    dev = getBackend(Val(args.device))
    # MPM procedure
    resetgridstatus_MLS_OS!(dev)(ndrange=grid.ni, grid)
    resetmpstatus_MLS_OS!(dev)(ndrange=mp.np, grid, mp)
    # F-bar based volumetric locking elimination approach
    if args.MVL == false
        aUpdatestatus_OS!(dev)(ndrange=mp.np, mp, ΔT)
    else
        aUpdatestatusvl1_OS!(dev)(ndrange=mp.np, grid, mp, ΔT)
        fastdiv!(dev)(ndrange=grid.ni, grid)
        aUpdatestatusvl2_OS!(dev)(ndrange=mp.np, grid, mp, ΔT)
    end
    # Update stress status
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
    # MPM procedure
    aP2G_MLS_OS!(dev)(ndrange=mp.np, grid, mp, ΔT)
    solvegrid_MLS_OS!(dev)(ndrange=grid.ni, grid, bc, G, ΔT)
    aG2P_MLS_OS!(dev)(ndrange=mp.np, grid, mp, attr, ΔT)
    return nothing
end