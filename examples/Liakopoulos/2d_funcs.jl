#==========================================================================================+
| Extension struct for the traction boundary                                               |
+==========================================================================================#
struct TractionBoundary{T<:AbstractArray} <: UserBoundaryExtra
    id::T
end

@user_struct TractionBoundary
#==========================================================================================#

@kernel function tP2G_TS!(
    grid   ::    DeviceGrid2D{T1, T2},
    mp     ::DeviceParticle2D{T1, T2},
    gravity::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        vol, ns, nl, σw = mp.Ω[ix], 1 - mp.n[ix], mp.S[ix] * mp.n[ix], mp.σw[ix]
        mps, mpw, mpm = mp.ms[ix], mp.mw[ix], mp.ms[ix] * ns + mp.mw[ix] * nl
        mppsx, mppsy, mppwx, mppwy = mp.ps[ix, 1], mp.ps[ix, 2], mp.pw[ix, 1], mp.pw[ix, 2]
        mpvsx, mpvsy, mpvwx, mpvwy = mp.vs[ix, 1], mp.vs[ix, 2], mp.vw[ix, 1], mp.vw[ix, 2]
        drag = (nl * mpw * T2(9.8)) / (mp.k[ix] * 4.41e-6)
        σxx, σyy, σxy = mp.σij[ix, 1], mp.σij[ix, 2], mp.σij[ix, 4]
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Nij = mp.Nij[ix, iy]
            if Nij ≠ T2(0.0)
                ∂Nx = mp.∂Nx[ix, iy]
                ∂Ny = mp.∂Ny[ix, iy]
                p2n = mp.p2n[ix, iy]
                # compute nodal mass for liquid & solid
                @KAatomic grid.mi[p2n] += Nij * mpw
                @KAatomic grid.mw[p2n] += Nij * mpw * nl
                @KAatomic grid.ms[p2n] += Nij * mps * ns
                # compute nodal momentum for liquid & solid
                @KAatomic grid.pw[p2n, 1] += Nij * mppwx * nl
                @KAatomic grid.pw[p2n, 2] += Nij * mppwy * nl
                @KAatomic grid.ps[p2n, 1] += Nij * mppsx * ns
                @KAatomic grid.ps[p2n, 2] += Nij * mppsy * ns
                # compute nodal total force for liquid & solid
                @KAatomic grid.fw[p2n, 1] += -vol * ∂Nx * σw
                @KAatomic grid.fw[p2n, 2] += -vol * ∂Ny * σw + 
                                              Nij * mpw * gravity
                @KAatomic grid.fs[p2n, 1] += -vol * (∂Nx * (σxx + σw) + ∂Ny * σxy)
                @KAatomic grid.fs[p2n, 2] += -vol * (∂Ny * (σyy + σw) + ∂Nx * σxy) +
                                              Nij * mpm * gravity
                # compute nodal drag force
                @KAatomic grid.fd[p2n, 1] += Nij * drag * (mpvwx - mpvsx)
                @KAatomic grid.fd[p2n, 2] += Nij * drag * (mpvwy - mpvsy)
            end
        end
    end
end

@kernel function tsolvegrid_TS!(
    grid::     DeviceGrid2D{T1, T2},
    bc  ::DeviceVBoundary2D{T1, T2},
    ΔT  ::T2,
    ζs  ::T2,
    ζw  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        mi_denom = grid.mi[ix] < eps(T2) ? T2(0.0) : inv(grid.mi[ix])
        mw_denom = grid.mw[ix] < eps(T2) ? T2(0.0) : inv(grid.mw[ix])
        # boundary condition
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vx_w_idx[ix] ≠ T1(0) ? grid.pw[ix, 1] = bc.vx_w_val[ix] : nothing
        bc.vy_w_idx[ix] ≠ T1(0) ? grid.pw[ix, 2] = bc.vy_w_val[ix] : nothing
        # compute nodal velocity
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        grid.vw[ix, 1] = grid.pw[ix, 1] * mw_denom
        grid.vw[ix, 2] = grid.pw[ix, 2] * mw_denom
        # compute damping force
        dampvw = -ζw * sqrt( grid.fw[ix, 1] * grid.fw[ix, 1]  + 
                             grid.fw[ix, 2] * grid.fw[ix, 2] )
        dampvs = -ζs * sqrt((grid.fs[ix, 1] - grid.fw[ix, 1]) * 
                            (grid.fs[ix, 1] - grid.fw[ix, 1]) +
                            (grid.fs[ix, 2] - grid.fw[ix, 2]) * 
                            (grid.fs[ix, 2] - grid.fw[ix, 2]))
        # compute node acceleration
        awx = mi_denom * (grid.fw[ix, 1] + dampvw * sign(grid.vw[ix, 1]) - grid.fd[ix, 1])
        awy = mi_denom * (grid.fw[ix, 2] + dampvw * sign(grid.vw[ix, 2]) - grid.fd[ix, 2])
        asx = ms_denom * (-grid.mw[ix] * awx + grid.fs[ix, 1] + 
            dampvw * sign(grid.vw[ix, 1]) + dampvs * sign(grid.vs[ix, 1]))
        asy = ms_denom * (-grid.mw[ix] * awy + grid.fs[ix, 2] + 
            dampvw * sign(grid.vw[ix, 2]) + dampvs * sign(grid.vs[ix, 2]))
        # update nodal temp velocity
        grid.vsT[ix, 1] = grid.vs[ix, 1] + asx * ΔT
        grid.vsT[ix, 2] = grid.vs[ix, 2] + asy * ΔT
        grid.vwT[ix, 1] = grid.vw[ix, 1] + awx * ΔT
        grid.vwT[ix, 2] = grid.vw[ix, 2] + awy * ΔT
        # apply boundary condition
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vx_w_idx[ix] ≠ T1(0) ? grid.vwT[ix, 1] = bc.vx_w_val[ix] : nothing
        bc.vy_w_idx[ix] ≠ T1(0) ? grid.vwT[ix, 2] = bc.vy_w_val[ix] : nothing
        # reset grid momentum
        grid.ps[ix, 1] = T2(0.0)
        grid.ps[ix, 2] = T2(0.0)
        grid.pw[ix, 1] = T2(0.0)
        grid.pw[ix, 2] = T2(0.0)
    end
end

@kernel function tdoublemapping1_TS!(
    grid    ::    DeviceGrid2D{T1, T2},
    mp      ::DeviceParticle2D{T1, T2},
    attr    ::  DeviceProperty{T1, T2},
    ΔT      ::T2,
    FLIP    ::T2,
    PIC     ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        ξsx = ξsy = ξwx = ξwy = vsx = vsy = vwx = vwy = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                ξsx += Ni *  grid.vsT[p2n, 1]
                ξsy += Ni *  grid.vsT[p2n, 2]
                ξwx += Ni *  grid.vwT[p2n, 1]
                ξwy += Ni *  grid.vwT[p2n, 2]
                vsx += Ni * (grid.vsT[p2n, 1] - grid.vs[p2n, 1])
                vsy += Ni * (grid.vsT[p2n, 2] - grid.vs[p2n, 2])
                vwx += Ni * (grid.vwT[p2n, 1] - grid.vw[p2n, 1])
                vwy += Ni * (grid.vwT[p2n, 2] - grid.vw[p2n, 2])
            end
        end
        # update particle ξition
        mp.ξ[ix, 1] += ΔT * ξsx
        mp.ξ[ix, 2] += ΔT * ξsy
        # update particle velocity
        mp.vs[ix, 1] = FLIP * (mp.vs[ix, 1] + vsx) + PIC * ξsx
        mp.vs[ix, 2] = FLIP * (mp.vs[ix, 2] + vsy) + PIC * ξsy
        mp.vw[ix, 1] = FLIP * (mp.vw[ix, 1] + vwx) + PIC * ξwx
        mp.vw[ix, 2] = FLIP * (mp.vw[ix, 2] + vwy) + PIC * ξwy
        # update particle momentum
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1]
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2]
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2]
    end
end

@kernel function tdoublemapping2_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # get particle porosity for liquid & solid
        nl, ns = mp.S[ix] * mp.n[ix], T2(1.0) - mp.n[ix]
        mppwx, mppwy, mppsx, mppsy = mp.pw[ix, 1], mp.pw[ix, 2], mp.ps[ix, 1], mp.ps[ix, 2]
        # update grid momentum for liquid & solid
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Nij = mp.Nij[ix, iy]
            if Nij ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                @KAatomic grid.pw[p2n, 1] += Nij * mppwx * nl
                @KAatomic grid.pw[p2n, 2] += Nij * mppwy * nl
                @KAatomic grid.ps[p2n, 1] += Nij * mppsx * ns
                @KAatomic grid.ps[p2n, 2] += Nij * mppsy * ns
            end
        end
    end
end

@kernel function tdoublemapping3_TS!(
    grid::     DeviceGrid2D{T1, T2},
    bc  ::DeviceVBoundary2D{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        mw_denom = grid.mw[ix] < eps(T2) ? T2(0.0) : inv(grid.mw[ix])
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        # compute nodal velocities for liquid & solid
        grid.vw[ix, 1] = grid.pw[ix, 1] * mw_denom
        grid.vw[ix, 2] = grid.pw[ix, 2] * mw_denom
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        # fixed Dirichlet nodes for liquid & solid
        bc.vx_w_idx[ix] ≠ T1(0) ? grid.vw[ix, 1] = bc.vx_w_val[ix] : nothing
        bc.vy_w_idx[ix] ≠ T1(0) ? grid.vw[ix, 2] = bc.vy_w_val[ix] : nothing
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 2] = bc.vy_s_val[ix] : nothing
        # compute nodal displacement for liquid & solid
        grid.Δuw[ix, 1] = grid.vw[ix, 1] * ΔT
        grid.Δuw[ix, 2] = grid.vw[ix, 2] * ΔT
        grid.Δus[ix, 1] = grid.vs[ix, 1] * ΔT
        grid.Δus[ix, 2] = grid.vs[ix, 2] * ΔT
    end
end

@kernel function tG2P_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
    attr::  DeviceProperty{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        ΔT_1 = inv(ΔT)
        dfs1 = dfs2 = dfs3 = dfs4 = dfw1 = dfw2 = dfw3 = dfw4 = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            if mp.Nij[ix, iy] ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                ∂Nx = mp.∂Nx[ix, iy]
                ∂Ny = mp.∂Ny[ix, iy]
                # compute solid incremental deformation gradient
                dfs1 += grid.Δus[p2n, 1] * ∂Nx
                dfs2 += grid.Δus[p2n, 1] * ∂Ny
                dfs3 += grid.Δus[p2n, 2] * ∂Nx
                dfs4 += grid.Δus[p2n, 2] * ∂Ny
                dfw1 += grid.Δuw[p2n, 1] * ∂Nx
                dfw2 += grid.Δuw[p2n, 1] * ∂Ny
                dfw3 += grid.Δuw[p2n, 2] * ∂Nx
                dfw4 += grid.Δuw[p2n, 2] * ∂Ny
            end
        end
        mp.ΔFs[ix, 1] = dfs1
        mp.ΔFs[ix, 2] = dfs2
        mp.ΔFs[ix, 3] = dfs3
        mp.ΔFs[ix, 4] = dfs4
        # strain rate (Second Invariant of Strain Rate Tensor)
        ϵvxx = dfs1 * ΔT_1
        ϵvyy = dfs4 * ΔT_1
        ϵvxy = T2(0.5) * (dfs2 + dfs3) * ΔT_1
        mp.ϵv[ix] = sqrt(ϵvxx * ϵvxx + ϵvyy * ϵvyy + T2(2.0) * ϵvxy * ϵvxy)
        # compute strain increment 
        mp.Δϵijs[ix, 1] = dfs1
        mp.Δϵijs[ix, 2] = dfs4
        mp.Δϵijs[ix, 4] = dfs2 + dfs3
        mp.Δϵijw[ix, 1] = dfw1
        mp.Δϵijw[ix, 2] = dfw4
        mp.Δϵijw[ix, 4] = dfw2 + dfw3
        # update strain tensor
        mp.ϵijs[ix, 1] += dfs1
        mp.ϵijs[ix, 2] += dfs4
        mp.ϵijs[ix, 4] += dfs2 + dfs3
        mp.ϵijw[ix, 1] += dfw1
        mp.ϵijw[ix, 2] += dfw4
        mp.ϵijw[ix, 4] += dfw2 + dfw3
        # deformation gradient matrix
        F1 = mp.F[ix, 1]; F2 = mp.F[ix, 2]; F3 = mp.F[ix, 3]; F4 = mp.F[ix, 4]      
        mp.F[ix, 1] = (dfs1 + T2(1.0)) * F1 + dfs2 * F3
        mp.F[ix, 2] = (dfs1 + T2(1.0)) * F2 + dfs2 * F4
        mp.F[ix, 3] = (dfs4 + T2(1.0)) * F3 + dfs3 * F1
        mp.F[ix, 4] = (dfs4 + T2(1.0)) * F4 + dfs3 * F2
        # update jacobian value and particle volume
        J         = mp.F[ix, 1] * mp.F[ix, 4] - mp.F[ix, 2] * mp.F[ix, 3]
        ΔJ        = J * mp.Ω0[ix] / mp.Ω[ix]
        mp.Ω[ix]  = J * mp.Ω0[ix]
        mp.ρs[ix] = mp.ρs0[ix] / J
        mp.ρw[ix] = mp.ρw0[ix] / J

        mp.n[ix] = clamp(T2(1.0) - (T2(1.0) - mp.n[ix]) / ΔJ, T2(0.0), T2(1.0))

        # update pore pressure
        suction = -mp.σw[ix]
        if suction ≥ T2(0.0)
            mp.σw[ix] += (attr.Kw[attr.nid[ix]] / mp.n[ix]) * (
                (T2(1.0) - mp.n[ix]) * (dfs1 + dfs4) + 
                           mp.n[ix]  * (dfw1 + dfw4))

            mp.S[ix] = T2(1 - 0.10152 * (suction / 9800)^2.4279)
            mp.k[ix] = T2(1 - 2.207 * (1 - mp.S[ix]))
            Cs
        else
            mp.σw[ix] += inv(Cs + mp.n[ix] * mp.S[ix]/attr.Kw[attr.nid[ix]]) * (
                (T2(1.0) - mp.n[ix]) * mp.S[ix] * (dfs1 + dfs4) + 
                           mp.n[ix]  * mp.S[ix] * (dfw1 + dfw4))
        end
        suction = -mp.σw[ix]
    end
end

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
    tP2G_TS!(dev)(ndrange=mp.np, grid, mp, G)

    # traction = T2(-1e3) / length(bc.ext.id)
    # grid.fs[bc.ext.id, 2] .+= traction

    tsolvegrid_TS!(dev)(ndrange=grid.ni, grid, bc, ΔT, args.ζs, args.ζw)
    tdoublemapping1_TS!(dev)(ndrange=mp.np, grid, mp, attr, ΔT, args.FLIP, args.PIC)
    tdoublemapping2_TS!(dev)(ndrange=mp.np, grid, mp)
    tdoublemapping3_TS!(dev)(ndrange=grid.ni, grid, bc, ΔT)
    tG2P_TS!(dev)(ndrange=mp.np, grid, mp, attr, ΔT)

    liE!(dev)(ndrange=mp.np, mp, attr)
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