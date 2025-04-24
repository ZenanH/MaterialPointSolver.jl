using KernelAbstractions

@kernel function tresetgridstatus_TS!(
    grid::DeviceGrid2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        if ix ≤ grid.nc
            grid.σm[ix] = T2(0.0)
            grid.σw[ix] = T2(0.0)
            grid.Ω[ix]  = T2(0.0)
        end
        grid.ms[ix]     = T2(0.0)
        grid.mi[ix]     = T2(0.0)
        grid.mw[ix]     = T2(0.0)
        grid.ps[ix, 1]  = T2(0.0)
        grid.ps[ix, 2]  = T2(0.0)
        grid.pw[ix, 1]  = T2(0.0)
        grid.pw[ix, 2]  = T2(0.0)
        grid.fs[ix, 1]  = T2(0.0)
        grid.fs[ix, 2]  = T2(0.0)
        grid.fw[ix, 1]  = T2(0.0)
        grid.fw[ix, 2]  = T2(0.0)
        grid.fd[ix, 1]  = T2(0.0)
        grid.fd[ix, 2]  = T2(0.0)
        #grid.ext.nl[ix] = T2(0.0)
    end
end

@kernel function tP2G_TS!(
    grid   ::    DeviceGrid2D{T1, T2},
    mp     ::DeviceParticle2D{T1, T2},
    gravity::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        S, vol = mp.S[ix], mp.Ω[ix]
        ns, nl = T2(1.0) - mp.n[ix], mp.n[ix] * S
        mps, mpw, mpm = mp.ms[ix], mp.mw[ix], mp.ms[ix] * ns + mp.mw[ix] * nl
        mppsx, mppsy, mppwx, mppwy = mp.ps[ix, 1], mp.ps[ix, 2], mp.pw[ix, 1], mp.pw[ix, 2]
        mpvsx, mpvsy, mpvwx, mpvwy = mp.vs[ix, 1], mp.vs[ix, 2], mp.vw[ix, 1], mp.vw[ix, 2]
        drag = (nl * mpw * T2(9.8)) / mp.k[ix]
        σxx, σyy, σxy, σw = mp.σij[ix, 1], mp.σij[ix, 2], mp.σij[ix, 4], mp.σw[ix]
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Nij = mp.Nij[ix, iy]
            if Nij ≠ T2(0.0)
                ∂Nx = mp.∂Nx[ix, iy]
                ∂Ny = mp.∂Ny[ix, iy]
                p2n = mp.p2n[ix, iy]
                # compute nodal mass for liquid & solid
                @KAatomic grid.mw[p2n] += Nij * mpw
                @KAatomic grid.mi[p2n] += Nij * mpw * nl
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
                @KAatomic grid.fs[p2n, 1] += -vol * (∂Nx * (σxx + S*σw) + ∂Ny * σxy)
                @KAatomic grid.fs[p2n, 2] += -vol * (∂Ny * (σyy + S*σw) + ∂Nx * σxy) +
                                              Nij * mpm * gravity
                # compute nodal drag force
                @KAatomic grid.fd[p2n, 1] += Nij * drag * (mpvwx - mpvsx)
                @KAatomic grid.fd[p2n, 2] += Nij * drag * (mpvwy - mpvsy)
                # compute nodel nl, e.i., liquid volume fraction
                #@KAatomic grid.ext.nl[p2n] += Nij * nl
            end
        end
    end
end

@kernel function tsolvegrid_TS!(
    grid::     DeviceGrid2D{T1, T2},
    attr::   DeviceProperty{T1, T2},
    bc  ::DeviceVBoundary2D{T1, T2},
    ΔT  ::T2,
    ζs  ::T2,
    ζw  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        mi_denom = grid.mi[ix] < eps(T2) ? T2(0.0) : inv(grid.mi[ix])
        mw_denom = grid.mw[ix] < eps(T2) ? T2(0.0) : inv(grid.mw[ix])
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        mt_denom = grid.ms[ix] + grid.mi[ix] < eps(T2) ? T2(0.0) : 
            inv((grid.ms[ix] + grid.mi[ix]) * ΔT)
        # boundary condition
        bc.vx_w_idx[ix] ≠ T1(0) ? grid.pw[ix, 1] = bc.vx_w_val[ix] : nothing
        bc.vy_w_idx[ix] ≠ T1(0) ? grid.pw[ix, 2] = bc.vy_w_val[ix] : nothing
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 2] = bc.vy_s_val[ix] : nothing
        # compute nodal velocity
        grid.vw[ix, 1] = grid.pw[ix, 1] * mi_denom
        grid.vw[ix, 2] = grid.pw[ix, 2] * mi_denom
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        # compute damping force
        dampvw = -ζw * sqrt( grid.fw[ix, 1] * grid.fw[ix, 1]  + 
                             grid.fw[ix, 2] * grid.fw[ix, 2] )
        dampvs = -ζs * sqrt((grid.fs[ix, 1] - grid.fw[ix, 1]) * 
                            (grid.fs[ix, 1] - grid.fw[ix, 1]) +
                            (grid.fs[ix, 2] - grid.fw[ix, 2]) * 
                            (grid.fs[ix, 2] - grid.fw[ix, 2]))
        # compute node acceleration
        awx = mw_denom * (grid.fw[ix, 1] + dampvw * sign(grid.vw[ix, 1]) - grid.fd[ix, 1])
        awy = mw_denom * (grid.fw[ix, 2] + dampvw * sign(grid.vw[ix, 2]) - grid.fd[ix, 2])
        asx = ms_denom * (-grid.mi[ix] * awx + grid.fs[ix, 1] + 
            dampvw * sign(grid.vw[ix, 1]) + dampvs * sign(grid.vs[ix, 1]))
        asy = ms_denom * (-grid.mi[ix] * awy + grid.fs[ix, 2] + 
            dampvw * sign(grid.vw[ix, 2]) + dampvs * sign(grid.vs[ix, 2]))
        # update nodal temp velocity
        grid.vwT[ix, 1] = grid.vw[ix, 1] + awx * ΔT
        grid.vwT[ix, 2] = grid.vw[ix, 2] + awy * ΔT
        grid.vsT[ix, 1] = grid.vs[ix, 1] + asx * ΔT
        grid.vsT[ix, 2] = grid.vs[ix, 2] + asy * ΔT
        # apply boundary condition
        bc.vx_w_idx[ix] ≠ T1(0) ? grid.vwT[ix, 1] = bc.vx_w_val[ix] : nothing
        bc.vy_w_idx[ix] ≠ T1(0) ? grid.vwT[ix, 2] = bc.vy_w_val[ix] : nothing
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 2] = bc.vy_s_val[ix] : nothing


        #==================================================================================+
        | Infiltration boundary for rainfall                                               |
        ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓=#
        # if ix in bc.ext.id
        #     # liquid volumetric fraction, infiltration velocity, and unit external normal 
        #     # vector
        #     nl, Ŵ, nx, ny = grid.ext.nl[ix], attr.ext.Ŵ, T2(0.0), T2(1.0)
        #     # net infiltration discharge
        #     Q_net_x = (nl * (grid.vwT[ix, 1] - grid.vsT[ix, 1]) - Ŵ) * nx
        #     Q_net_y = (nl * (grid.vwT[ix, 2] - grid.vsT[ix, 2]) - Ŵ) * ny
        #     # velocity correction on x direction
        #     # if Q_net_x ≥ T2(0.0)
        #     #     awx = (grid.fs[ix, 1] * ΔT - grid.ms[ix] * (
        #     #         grid.vw[ix, 1] - grid.vs[ix, 1] - grid.vwT[ix, 1] + grid.vsT[ix, 1])
        #     #         ) / ((grid.mi[ix] + grid.ms[ix]) * ΔT)
        #     #     asx = (grid.fs[ix, 1] - grid.mi[ix] * awx) / grid.ms[ix]
        #     #     # update nodal temp velocity
        #     #     grid.vwT[ix, 1] = grid.vw[ix, 1] + awx * ΔT
        #     #     #grid.vsT[ix, 1] = grid.vs[ix, 1] + asx * ΔT
        #     # end
        #     # velocity correction on y direction
        #     if Q_net_y ≤ T2(0.0)
        #         # awy = (grid.fs[ix, 2] * ΔT - grid.ms[ix] * (
        #         #     grid.vw[ix, 2] - grid.vs[ix, 2] - grid.vwT[ix, 2] + grid.vsT[ix, 2])
        #         #     ) * mt_denom
        #         # asy = (grid.fs[ix, 2] - grid.mi[ix] * awy) * ms_denom
        #         # # update nodal temp velocity
        #         # grid.vwT[ix, 2] = grid.vw[ix, 2] + awy * ΔT
        #         # grid.vsT[ix, 2] = grid.vs[ix, 2] + asy * ΔT
        #         #grid.vwT[ix, 2] = Ŵ * inv(nl) + grid.vsT[ix, 2]
        #     end
        #     # apply boundary condition
        #     bc.vx_w_idx[ix] ≠ T1(0) ? grid.vwT[ix, 1] = bc.vx_w_val[ix] : nothing
        #     bc.vy_w_idx[ix] ≠ T1(0) ? grid.vwT[ix, 2] = bc.vy_w_val[ix] : nothing
        #     bc.vx_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 1] = bc.vx_s_val[ix] : nothing
        #     bc.vy_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 2] = bc.vy_s_val[ix] : nothing
        # end
        #↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑#

        # reset grid momentum
        grid.pw[ix, 1] = T2(0.0)
        grid.pw[ix, 2] = T2(0.0)
        grid.ps[ix, 1] = T2(0.0)
        grid.ps[ix, 2] = T2(0.0)
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
        # prepare variables to accumulate the nodal values
        ξsx = ξsy = ξwx = ξwy = vsx = vsy = vwx = vwy = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                p2n  = mp.p2n[ix, iy]
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
        # mp.ξ[ix, 1] += ΔT * ξsx
        # mp.ξ[ix, 2] += ΔT * ξsy
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
        nl, ns = mp.n[ix]*mp.S[ix], T2(1.0) - mp.n[ix]
        # update grid momentum for liquid & solid
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Nij = mp.Nij[ix, iy]
            if Nij ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                @KAatomic grid.pw[p2n, 1] += mp.pw[ix, 1] * Nij * nl
                @KAatomic grid.pw[p2n, 2] += mp.pw[ix, 2] * Nij * nl
                @KAatomic grid.ps[p2n, 1] += mp.ps[ix, 1] * Nij * ns
                @KAatomic grid.ps[p2n, 2] += mp.ps[ix, 2] * Nij * ns
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
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        mi_denom = grid.mi[ix] < eps(T2) ? T2(0.0) : inv(grid.mi[ix])
        # compute nodal velocities for liquid & solid
        grid.vw[ix, 1] = grid.pw[ix, 1] * mi_denom
        grid.vw[ix, 2] = grid.pw[ix, 2] * mi_denom
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
        # get volumetric fraction parameters
        nl, ni = mp.n[ix] * mp.S[ix], mp.S[ix] * (T2(1.0) - mp.n[ix])
        # update pore pressure
        σw = mp.σw[ix] < 0 ? -mp.σw[ix] : mp.σw[ix]
        Pc = clamp(σw, T2(0.0), T2(18154.0))
    
        ∂S∂s = T2(-5.0241518982722e-11)*Pc^T2(1.4279)
        mp.σw[ix] += inv(mp.n[ix] * -∂S∂s + nl * inv(attr.Kw[attr.nid[ix]])) * (
            nl * (dfw1 + dfw4) + ni * (dfs1 + dfs4))
        # update saturation by SWCC
        mp.S[ix] = T2(1.0) - T2(0.10152) * (Pc * T2(0.000102))^T2(2.4279)
        # update hydraulic conductivity by HCC
        mp.k[ix] = (T2(1.0) - T2(2.207) * (T2(1.0) - mp.S[ix])) * T2(4.41e-6)
        # update porosity
        # mp.n[ix] = clamp(T2(1.0) - (T2(1.0) - mp.n[ix]) / ΔJ, T2(0.0), T2(1.0))
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
    # MPM procedure
    tresetgridstatus_TS!(dev)(ndrange=grid.ni, grid)
    resetmpstatus_TS!(dev)(ndrange=mp.np, grid, mp, Val(args.basis))
    tP2G_TS!(dev)(ndrange=mp.np, grid, mp, G)
    tsolvegrid_TS!(dev)(ndrange=grid.ni, grid, attr, bc, ΔT, args.ζs, args.ζw)

    #grid.vwT[bc.ext.id, 2] .= attr.ext.Ŵy

    tdoublemapping1_TS!(dev)(ndrange=mp.np, grid, mp, attr, ΔT, args.FLIP, args.PIC)
    tdoublemapping2_TS!(dev)(ndrange=mp.np, grid, mp)
    tdoublemapping3_TS!(dev)(ndrange=grid.ni, grid, bc, ΔT)
    tG2P_TS!(dev)(ndrange=mp.np, grid, mp, attr, ΔT)
    # update stress status
    liE!(dev)(ndrange=mp.np, mp, attr)
    # cell-averaged volumetric locking elimination
    if args.MVL == true
        vollock1_TS!(dev)(ndrange=mp.np, grid, mp)
        vollock2_TS!(dev)(ndrange=mp.np, grid, mp)
    end
    return nothing
end