using KernelAbstractions

@kernel function tresetgridstatus_TS!(
    grid::DeviceGrid2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        grid.σm[ix]    = T2(0.0)
        grid.σw[ix]    = T2(0.0)
        grid.Ω[ix]     = T2(0.0)
        grid.ms[ix]    = T2(0.0)
        grid.mi[ix]    = T2(0.0)
        grid.mw[ix]    = T2(0.0)
        grid.ps[ix, 1] = T2(0.0)
        grid.ps[ix, 2] = T2(0.0)
        grid.pw[ix, 1] = T2(0.0)
        grid.pw[ix, 2] = T2(0.0)
        grid.fs[ix, 1] = T2(0.0)
        grid.fs[ix, 2] = T2(0.0)
        grid.fw[ix, 1] = T2(0.0)
        grid.fw[ix, 2] = T2(0.0)
        grid.fd[ix, 1] = T2(0.0)
        grid.fd[ix, 2] = T2(0.0)
    end
end

@kernel function G2Pvl1_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        dfs1 = dfs2 = dfs3 = dfs4 = T2(0.0)
        dfw1 = dfw2 = dfw3 = dfw4 = T2(0.0)
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
        dfs1 += T2(1.0); dfs4 += T2(1.0)
        dfw1 += T2(1.0); dfw4 += T2(1.0)
        mp.ΔFs[ix, 1] = dfs1; mp.ΔFs[ix, 2] = dfs2 
        mp.ΔFs[ix, 3] = dfs3; mp.ΔFs[ix, 4] = dfs4
        mp.ΔFw[ix, 1] = dfw1; mp.ΔFw[ix, 2] = dfw2
        mp.ΔFw[ix, 3] = dfw3; mp.ΔFw[ix, 4] = dfw4
        # compute ΔJₚ in the current time step
        ΔJs = dfs1 * dfs4 - dfs2 * dfs3
        ΔJw = dfw1 * dfw4 - dfw2 * dfw3 
        # map this value from particle to grid cell
        vol = mp.Ω[ix]
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            NiV = mp.Nij[ix, iy] * vol
            if NiV ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                @KAatomic grid.σm[p2n] += NiV * ΔJs
                @KAatomic grid.σw[p2n] += NiV * ΔJw
                @KAatomic grid.Ω[p2n] += NiV
            end
        end
    end
end

@kernel function fastdiv_TS!(grid::DeviceGrid{T1, T2}) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        Ω = grid.Ω[ix] == T2(0.0) ? T2(0.0) : inv(grid.Ω[ix])
        grid.σm[ix] *= Ω
        grid.σw[ix] *= Ω
    end
end

@kernel function G2Pvl2_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
    attr::  DeviceProperty{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    ΔT_1 = inv(ΔT)
    if ix ≤ mp.np
        cos, cow = T2(0.0), T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Nij = mp.Nij[ix, iy]
            if Nij ≠ T2(0.0)
                p2n  = mp.p2n[ix, iy]
                cos += Nij * grid.σm[p2n]
                cow += Nij * grid.σw[p2n]
            end
        end
        Jrs = mp.ΔFs[ix, 1] * mp.ΔFs[ix, 4] - mp.ΔFs[ix, 2] * mp.ΔFs[ix, 3]
        Jrw = mp.ΔFw[ix, 1] * mp.ΔFw[ix, 4] - mp.ΔFw[ix, 2] * mp.ΔFw[ix, 3]
        Jcs = (cos / Jrs) ^ T2(0.5)
        Jcw = (cow / Jrw) ^ T2(0.5)
        # update deformation gradient matrix
        mp.ΔFs[ix, 1] *= Jcs; mp.ΔFs[ix, 2] *= Jcs
        mp.ΔFs[ix, 3] *= Jcs; mp.ΔFs[ix, 4] *= Jcs
        mp.ΔFw[ix, 1] *= Jcw; mp.ΔFw[ix, 2] *= Jcw
        mp.ΔFw[ix, 3] *= Jcw; mp.ΔFw[ix, 4] *= Jcw
        mp.ΔFs[ix, 1] -= T2(1.0); mp.ΔFs[ix, 4] -= T2(1.0)
        mp.ΔFw[ix, 1] -= T2(1.0); mp.ΔFw[ix, 4] -= T2(1.0)
        dfs1 = mp.ΔFs[ix, 1]; dfs2 = mp.ΔFs[ix, 2]
        dfs3 = mp.ΔFs[ix, 3]; dfs4 = mp.ΔFs[ix, 4]
        dfw1 = mp.ΔFw[ix, 1]; dfw2 = mp.ΔFw[ix, 2]
        dfw3 = mp.ΔFw[ix, 3]; dfw4 = mp.ΔFw[ix, 4]
        # strain rate (Second Invariant of Strain Rate Tensor)
        dϵxx = dfs1 * ΔT_1
        dϵyy = dfs4 * ΔT_1
        dϵxy = T2(0.5) * (dfs2 + dfs3) * ΔT_1
        mp.ϵv[ix] = sqrt(dϵxx * dϵxx + dϵyy * dϵyy + T2(2.0) * dϵxy * dϵxy)
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
        # update jacobian value and particle Ωume
        J         = mp.F[ix, 1] * mp.F[ix, 4] - mp.F[ix, 2] * mp.F[ix, 3]
        ΔJ        = J * mp.Ω0[ix] / mp.Ω[ix]
        mp.Ω[ix]  = J * mp.Ω0[ix]
        mp.ρs[ix] = mp.ρs0[ix] / J
        mp.ρw[ix] = mp.ρw0[ix] / J
        # update pore pressure and n
        mp.σw[ix] += (attr.Kw[attr.nid[ix]] / mp.n[ix]) * (
            (T2(1.0) - mp.n[ix]) * (dfs1 + dfs4) + 
                       mp.n[ix]  * (dfw1 + dfw4))
        mp.n[ix] = clamp(T2(1.0) - (T2(1.0) - mp.n[ix]) / ΔJ, T2(0.0), T2(1.0))
    end
end

function tprocedure!(
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
    G = Ti < args.Te ? args.gravity / args.Te * Ti : args.gravity
    dev = getBackend(Val(args.device))
    # MPM procedure
    tresetgridstatus_TS!(dev)(ndrange=grid.ni, grid)
    resetmpstatus_TS!(dev)(ndrange=mp.np, grid, mp, Val(args.basis))
    P2G_TS!(dev)(ndrange=mp.np, grid, mp, attr, G)
    solvegrid_TS!(dev)(ndrange=grid.ni, grid, bc, ΔT, args.ζs, args.ζw)
    doublemapping1_TS!(dev)(ndrange=mp.np, grid, mp, attr, ΔT, args.FLIP, args.PIC)
    doublemapping2_TS!(dev)(ndrange=mp.np, grid, mp)
    doublemapping3_TS!(dev)(ndrange=grid.ni, grid, bc, ΔT)
    # F-bar based volumetric locking elimination approach
    G2Pvl1_TS!(dev)(ndrange=mp.np, grid, mp)
    fastdiv_TS!(dev)(ndrange=grid.ni, grid)
    G2Pvl2_TS!(dev)(ndrange=mp.np, grid, mp, attr, ΔT)
    # update stress status
    liE!(dev)(ndrange=mp.np, mp, attr)
    return nothing
end