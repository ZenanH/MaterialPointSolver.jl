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
    P2G_TS!(dev)(ndrange=mp.np, grid, mp, G)

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

# cv = init_k / (9800 * (1 / init_Es + init_n / init_Kw))
# t = 0.1 / cv
# helper functions =====================================================================
function terzaghi(p0, Tv)
    num = 100
    H = 1
    Z = range(0, 1, length=num)
    data = zeros(num, 2)
    data[:, 2] .= Z
    @inbounds for i in 1:num
        p = 0.0
        for m in 1:2:1e4
            p += 4*p0/π*(1/m)*sin((m*π*data[i, 2])/(2*H))*exp((-m^2)*((π/2)^2)*Tv)
        end
        data[num+1-i, 1] = p/p0
    end
    return data
end

function consolidation()
    num = 1000
    dat = zeros(num, 2)    
    dat[:, 1] .= collect(range(0, 10, length=num))
    @inbounds for i in 1:num
        tmp = 0.0
        for m in 1:2:1e4
            tmp += (8/π^2)*(1/m^2)*exp(-(m*π/2)^2*dat[i, 1]) 
        end
        dat[i, 2] = 1-tmp
    end
    return dat
end