#==========================================================================================+
| Extension struct for the traction boundary                                               |
+==========================================================================================#
struct TractionBoundary{T<:AbstractArray} <: UserBoundaryExtra
    id::T
end

@user_struct TractionBoundary
#==========================================================================================#

function Tprocedure!(
    args::     DeviceArgs{T1, T2}, 
    grid::     DeviceGrid{T1, T2}, 
    mp  :: DeviceParticle{T1, T2}, 
    attr:: DeviceProperty{T1, T2},
    bc  ::DeviceVBoundary{T1, T2},
    ΔT  ::T2,
    Ti  ::T2,
        ::Val{:TS},
        ::Val{:MUSL}
) where {T1, T2}
    Ti < args.Te ? G = args.gravity / args.Te * Ti : G = args.gravity
    dev = getBackend(Val(args.device))

    resetgridstatus_TS!(dev)(ndrange=grid.ni, grid)
    resetmpstatus_TS!(dev)(ndrange=mp.np, grid, mp, Val(args.basis))
    P2G_TS!(dev)(ndrange=mp.np, grid, mp, attr, G)

    traction = T2(-1e3) / length(bc.ext.id)
    grid.fs[bc.ext.id, 2] .+= traction

    solvegrid_TS!(dev)(ndrange=grid.ni, grid, bc, ΔT, args.ζs, args.ζw)
    doublemapping1_TS!(dev)(ndrange=mp.np, grid, mp, attr, ΔT, args.FLIP, args.PIC)
    doublemapping2_TS!(dev)(ndrange=mp.np, grid, mp)
    doublemapping3_TS!(dev)(ndrange=grid.ni, grid, bc, ΔT)
    G2P_TS!(dev)(ndrange=mp.np, grid, mp, attr, ΔT)

    liE!(dev)(ndrange=mp.np, mp, attr)
    # volumetric locking elimination approach
    if args.MVL == true
        vollock1_TS!(dev)(ndrange=mp.np, grid, mp)
        vollock2_TS!(dev)(ndrange=mp.np, grid, mp)
    end
    return nothing
end