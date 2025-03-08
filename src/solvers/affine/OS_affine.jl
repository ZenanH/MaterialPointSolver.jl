#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : OS_affine.jl                                                               |
|  Description: Implementation of affine MPM                                               |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Functions  : 01. aUpdatestatus_OS!    [2D & 3D]                                         |
|               02. aP2G_OS!             [2D & 3D]                                         |
|               03. aG2P_OS!             [2D & 3D]                                         |
|               04. aUpdatestatusvl1_OS! [2D & 3D]                                         |
|               05. aUpdatestatusvl2_OS! [2D & 3D]                                         |
+==========================================================================================#

export aUpdatestatus_OS!
export aP2G_OS!
export aG2P_OS!
export aUpdatestatusvl1_OS!
export aUpdatestatusvl2_OS!

"""
    aUpdatestatus_OS!(mp::DeviceParticle2D{T1, T2}, ΔT::T2)

Description:
---
Update particle status: domain, stress, density ... (affine)
"""
@kernel inbounds = true function aUpdatestatus_OS!(
    mp     ::DeviceParticle2D{T1, T2},
    ΔT     ::T2
) where {T1, T2}
    ix = @index(Global)
    ΔT_1 = inv(ΔT)
    if ix ≤ mp.np
        # get values from affine matrix
        dF1, dF2 = mp.aC[ix, 1] * ΔT, mp.aC[ix, 2] * ΔT
        dF3, dF4 = mp.aC[ix, 3] * ΔT, mp.aC[ix, 4] * ΔT
        # update deformation gradient incremental matrix
        mp.ΔFs[ix, 1], mp.ΔFs[ix, 2] = dF1, dF2
        mp.ΔFs[ix, 3], mp.ΔFs[ix, 4] = dF3, dF4
        # strain rate (Second Invariant of Strain Rate Tensor)
        dϵxx = dF1 * ΔT_1
        dϵyy = dF4 * ΔT_1
        dϵxy = T2(0.5) * (dF2 + dF3) * ΔT_1
        mp.ϵv[ix] = sqrt(dϵxx * dϵxx + dϵyy * dϵyy + T2(2.0) * dϵxy * dϵxy)
        # update strain increment
        mp.Δϵijs[ix, 1] = dF1
        mp.Δϵijs[ix, 2] = dF4
        mp.Δϵijs[ix, 4] = dF2 + dF3
        # update strain tensor
        mp.ϵijs[ix, 1] += dF1
        mp.ϵijs[ix, 2] += dF4
        mp.ϵijs[ix, 4] += dF2 + dF3
        # deformation gradient matrix
        F1, F2, F3, F4 = mp.F[ix, 1], mp.F[ix, 2], mp.F[ix, 3], mp.F[ix, 4]  
        mp.F[ix, 1] = (dF1 + T2(1.0)) * F1 + dF2 * F3
        mp.F[ix, 2] = (dF1 + T2(1.0)) * F2 + dF2 * F4
        mp.F[ix, 3] = (dF4 + T2(1.0)) * F3 + dF3 * F1
        mp.F[ix, 4] = (dF4 + T2(1.0)) * F4 + dF3 * F2
        # update jacobian value and particle volume
        J = mp.F[ix, 1] * mp.F[ix, 4] - mp.F[ix, 2] * mp.F[ix, 3]
        mp.Ω[ix]  = J * mp.Ω0[ix]
        mp.ρs[ix] = mp.ρs0[ix] / J
    end
end

@kernel inbounds = true function aUpdatestatus_OS!(
    mp     ::DeviceParticle3D{T1, T2},
    ΔT     ::T2
) where {T1, T2}
    ix = @index(Global)
    ΔT_1 = inv(ΔT)
    if ix ≤ mp.np
        # get values from affine matrix
        dF1, dF2, dF3 = mp.aC[ix, 1] * ΔT, mp.aC[ix, 2] * ΔT, mp.aC[ix, 3] * ΔT
        dF4, dF5, dF6 = mp.aC[ix, 4] * ΔT, mp.aC[ix, 5] * ΔT, mp.aC[ix, 6] * ΔT
        dF7, dF8, dF9 = mp.aC[ix, 7] * ΔT, mp.aC[ix, 8] * ΔT, mp.aC[ix, 9] * ΔT
        # update deformation gradient incremental matrix
        mp.ΔFs[ix, 1], mp.ΔFs[ix, 2], mp.ΔFs[ix, 3] = dF1, dF2, dF3
        mp.ΔFs[ix, 4], mp.ΔFs[ix, 5], mp.ΔFs[ix, 6] = dF4, dF5, dF6
        mp.ΔFs[ix, 7], mp.ΔFs[ix, 8], mp.ΔFs[ix, 9] = dF7, dF8, dF9
        # strain rate (Second Invariant of Strain Rate Tensor)
        dϵxx = dF1 * ΔT_1
        dϵyy = dF5 * ΔT_1
        dϵzz = dF9 * ΔT_1
        dϵxy = T2(0.5) * (dF2 + dF4) * ΔT_1
        dϵyz = T2(0.5) * (dF6 + dF8) * ΔT_1
        dϵxz = T2(0.5) * (dF3 + dF7) * ΔT_1
        mp.ϵv[ix] = sqrt(dϵxx * dϵxx + dϵyy * dϵyy + dϵzz * dϵzz + 
            T2(2.0) * (dϵxy * dϵxy + dϵyz * dϵyz + dϵxz * dϵxz))
        # update strain increment
        mp.Δϵijs[ix, 1] = dF1
        mp.Δϵijs[ix, 2] = dF5
        mp.Δϵijs[ix, 3] = dF9
        mp.Δϵijs[ix, 4] = dF2 + dF4
        mp.Δϵijs[ix, 5] = dF6 + dF8
        mp.Δϵijs[ix, 6] = dF3 + dF7
        # update strain tensor
        mp.ϵijs[ix, 1] += dF1
        mp.ϵijs[ix, 2] += dF5
        mp.ϵijs[ix, 3] += dF9
        mp.ϵijs[ix, 4] += dF2 + dF4
        mp.ϵijs[ix, 5] += dF6 + dF8
        mp.ϵijs[ix, 6] += dF3 + dF7
        # deformation gradient matrix
        F1, F2, F3 = mp.F[ix, 1], mp.F[ix, 2], mp.F[ix, 3]
        F4, F5, F6 = mp.F[ix, 4], mp.F[ix, 5], mp.F[ix, 6]
        F7, F8, F9 = mp.F[ix, 7], mp.F[ix, 8], mp.F[ix, 9]
        mp.F[ix, 1] = (dF1 + T2(1.0)) * F1 + dF2 * F4 + dF3 * F7
        mp.F[ix, 2] = (dF1 + T2(1.0)) * F2 + dF2 * F5 + dF3 * F8
        mp.F[ix, 3] = (dF1 + T2(1.0)) * F3 + dF2 * F6 + dF3 * F9
        mp.F[ix, 4] = (dF5 + T2(1.0)) * F4 + dF4 * F1 + dF6 * F7
        mp.F[ix, 5] = (dF5 + T2(1.0)) * F5 + dF4 * F2 + dF6 * F8
        mp.F[ix, 6] = (dF5 + T2(1.0)) * F6 + dF4 * F3 + dF6 * F9
        mp.F[ix, 7] = (dF9 + T2(1.0)) * F7 + dF8 * F4 + dF7 * F1
        mp.F[ix, 8] = (dF9 + T2(1.0)) * F8 + dF8 * F5 + dF7 * F2
        mp.F[ix, 9] = (dF9 + T2(1.0)) * F9 + dF8 * F6 + dF7 * F3
        # update jacobian value and particle volume
        J = mp.F[ix, 1] * mp.F[ix, 5] * mp.F[ix, 9] + 
            mp.F[ix, 2] * mp.F[ix, 6] * mp.F[ix, 7] +
            mp.F[ix, 3] * mp.F[ix, 4] * mp.F[ix, 8] - 
            mp.F[ix, 7] * mp.F[ix, 5] * mp.F[ix, 3] -
            mp.F[ix, 8] * mp.F[ix, 6] * mp.F[ix, 1] - 
            mp.F[ix, 9] * mp.F[ix, 4] * mp.F[ix, 2] 
        mp.Ω[ix]  = J * mp.Ω0[ix]
        mp.ρs[ix] = mp.ρs0[ix] / J
    end
end

"""
    aP2G_OS!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2}, gravity::T2)

Description:
---
aP2G procedure for scattering the mass, momentum, and forces from particles to grid (affine).
"""
@kernel inbounds = true function aP2G_OS!(
    grid   ::    DeviceGrid2D{T1, T2},
    mp     ::DeviceParticle2D{T1, T2},
    gravity::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        vol = mp.Ω[ix]
        mps = mp.ms[ix]
        mppsx, mppsy = mp.ps[ix, 1], mp.ps[ix, 2]
        σxx, σyy, σxy = mp.σij[ix, 1], mp.σij[ix, 2], mp.σij[ix, 4]
        mξx, mξy = mp.ξ[ix, 1], mp.ξ[ix, 2]
        aC1, aC2, aC3, aC4 = mp.aC[ix, 1], mp.aC[ix, 2], mp.aC[ix, 3], mp.aC[ix, 4]
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                ∂Nx = mp.∂Nx[ix, iy]
                ∂Ny = mp.∂Ny[ix, iy]
                p2n = mp.p2n[ix, iy]
                dx, dy = grid.ξ[p2n, 1] - mξx, grid.ξ[p2n, 2] - mξy
                # compute nodal mass
                @KAatomic grid.ms[p2n] += mps * Ni
                # compute nodal momentum
                @KAatomic grid.ps[p2n, 1] += Ni * (mppsx + mps * (aC1 * dx + aC2 * dy))
                @KAatomic grid.ps[p2n, 2] += Ni * (mppsy + mps * (aC3 * dx + aC4 * dy))
                # compute nodal total force for solid
                @KAatomic grid.fs[p2n, 1] += -vol * (∂Nx * σxx + ∂Ny * σxy)
                @KAatomic grid.fs[p2n, 2] += -vol * (∂Ny * σyy + ∂Nx * σxy) + 
                    Ni * mps * gravity
            end
        end
    end
end

@kernel inbounds = true function aP2G_OS!(
    grid   ::    DeviceGrid3D{T1, T2},
    mp     ::DeviceParticle3D{T1, T2},
    gravity::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        mps, mx, my, mz, vol = mp.ms[ix], mp.ξ[ix, 1], mp.ξ[ix, 2], mp.ξ[ix, 3], mp.Ω[ix]
        mppx, mppy, mppz = mp.ps[ix, 1], mp.ps[ix, 2], mp.ps[ix, 3]
        aC1, aC2, aC3 = mp.aC[ix, 1], mp.aC[ix, 2], mp.aC[ix, 3]
        aC4, aC5, aC6 = mp.aC[ix, 4], mp.aC[ix, 5], mp.aC[ix, 6]
        aC7, aC8, aC9 = mp.aC[ix, 7], mp.aC[ix, 8], mp.aC[ix, 9]
        σxx, σyy, σzz = mp.σij[ix, 1], mp.σij[ix, 2], mp.σij[ix, 3]
        σxy, σyz, σzx = mp.σij[ix, 4], mp.σij[ix, 5], mp.σij[ix, 6]
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                ∂Nx = mp.∂Nx[ix, iy]
                ∂Ny = mp.∂Ny[ix, iy]
                ∂Nz = mp.∂Nz[ix, iy]
                p2n = mp.p2n[ix, iy]
                dx, dy, dz = grid.ξ[p2n, 1] - mx, grid.ξ[p2n, 2] - my, grid.ξ[p2n, 3] - mz
                # compute nodal mass
                @KAatomic grid.ms[p2n] += Ni * mps
                # compute nodal momentum
                @KAatomic grid.ps[p2n, 1] += Ni * (mppx + mps * 
                    (aC1 * dx + aC2 * dy + aC3 * dz))
                @KAatomic grid.ps[p2n, 2] += Ni * (mppy + mps *
                    (aC4 * dx + aC5 * dy + aC6 * dz))
                @KAatomic grid.ps[p2n, 3] += Ni * (mppz + mps *
                    (aC7 * dx + aC8 * dy + aC9 * dz))
                # compute nodal total force for solid
                @KAatomic grid.fs[p2n, 1] += -vol * (∂Nx * σxx + ∂Ny * σxy + ∂Nz * σzx)
                @KAatomic grid.fs[p2n, 2] += -vol * (∂Ny * σyy + ∂Nx * σxy + ∂Nz * σyz)
                @KAatomic grid.fs[p2n, 3] += -vol * (∂Nz * σzz + ∂Nx * σzx + ∂Ny * σyz) + 
                    Ni * mps * gravity
            end
        end
    end
end

"""
    aG2P_OS!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2}, 
        attr::DeviceProperty{T1, T2}, ΔT::T2)

Description:
---
Mapping results from grid to particles. (affine)
"""
@kernel inbounds = true function aG2P_OS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
    attr::  DeviceProperty{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        vx = vy = T2(0.0)
        B11 = B12 = B21 = B22 = T2(0.0)
        D11 = D12 = D22 = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            mx, my = mp.ξ[ix, 1], mp.ξ[ix, 2]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                gx, gy = grid.ξ[p2n, 1], grid.ξ[p2n, 2]
                vx += Ni * grid.vsT[p2n, 1]
                vy += Ni * grid.vsT[p2n, 2]
                # compute B matrix
                B11 += Ni * grid.vsT[p2n, 1] * (gx - mx)
                B12 += Ni * grid.vsT[p2n, 1] * (gy - my)
                B21 += Ni * grid.vsT[p2n, 2] * (gx - mx)
                B22 += Ni * grid.vsT[p2n, 2] * (gy - my)
                # compute D matrix
                D11 += Ni * (gx - mx) * (gx - mx)
                D12 += Ni * (gx - mx) * (gy - my)
                D22 += Ni * (gy - my) * (gy - my)
            end
        end
        # update particle velocity
        mp.vs[ix, 1] = vx
        mp.vs[ix, 2] = vy
        # update particle position
        mp.ξ[ix, 1] += ΔT * vx
        mp.ξ[ix, 2] += ΔT * vy
        # update affine matrix C
        D_del = D11 * D22 - D12 * D12
        D_del = D_del == T2(0.0) ? T2(1.0) : inv(D_del)
        d_11 =  D_del * D22
        d_12 = -D_del * D12
        d_22 =  D_del * D11
        mp.aC[ix, 1] = B11 * d_11 + B12 * d_12
        mp.aC[ix, 2] = B11 * d_12 + B12 * d_22
        mp.aC[ix, 3] = B21 * d_11 + B22 * d_12
        mp.aC[ix, 4] = B21 * d_12 + B22 * d_22
        # update CFL conditions
        nid = attr.nid[ix]
        Ks  = attr.Ks[nid]
        Gs  = attr.Gs[nid]
        cdil = sqrt((Ks + T2(1.333333) * Gs) / mp.ρs[ix]) # 4/3 ≈ 1.333333
        mp.cfl[ix] = min(grid.dx / (cdil + abs(mp.vs[ix, 1])), 
                         grid.dy / (cdil + abs(mp.vs[ix, 2]))) 
    end
end

@kernel inbounds = true function aG2P_OS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
    attr::  DeviceProperty{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        vx = vy = vz = T2(0.0)
        B11 = B12 = B13 = B21 = B22 = B23 = B31 = B32 = B33 = T2(0.0)
        D11 = D12 = D13 = D22 = D23 = D33 = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            mx, my, mz = mp.ξ[ix, 1], mp.ξ[ix, 2], mp.ξ[ix, 3]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                gx, gy, gz = grid.ξ[p2n, 1], grid.ξ[p2n, 2], grid.ξ[p2n, 3]
                vx += Ni * grid.vsT[p2n, 1]
                vy += Ni * grid.vsT[p2n, 2]
                vz += Ni * grid.vsT[p2n, 3]
                # compute B matrix
                B11 += Ni * grid.vsT[p2n, 1] * (gx - mx)
                B12 += Ni * grid.vsT[p2n, 1] * (gy - my)
                B13 += Ni * grid.vsT[p2n, 1] * (gz - mz)
                B21 += Ni * grid.vsT[p2n, 2] * (gx - mx)
                B22 += Ni * grid.vsT[p2n, 2] * (gy - my)
                B23 += Ni * grid.vsT[p2n, 2] * (gz - mz)
                B31 += Ni * grid.vsT[p2n, 3] * (gx - mx)
                B32 += Ni * grid.vsT[p2n, 3] * (gy - my)
                B33 += Ni * grid.vsT[p2n, 3] * (gz - mz)
                # compute D matrix
                D11 += Ni * (gx - mx) * (gx - mx)
                D12 += Ni * (gx - mx) * (gy - my)
                D13 += Ni * (gx - mx) * (gz - mz)
                D22 += Ni * (gy - my) * (gy - my)
                D23 += Ni * (gy - my) * (gz - mz)
                D33 += Ni * (gz - mz) * (gz - mz)
            end
        end
        # update particle velocity
        mp.vs[ix, 1] = vx
        mp.vs[ix, 2] = vy
        mp.vs[ix, 3] = vz
        # update particle position
        mp.ξ[ix, 1] += ΔT * vx
        mp.ξ[ix, 2] += ΔT * vy
        mp.ξ[ix, 3] += ΔT * vz
        # update affine matrix C
        D_del = D11 * D22 * D33 + D12 * D23 * D13 + D13 * D12 * D23 - 
                D13 * D22 * D13 - D12 * D12 * D33 - D11 * D23 * D23
        D_del = D_del == T2(0.0) ? T2(1.0) : inv(D_del)
        d_11 = D_del * (D22 * D33 - D23 * D23)
        d_12 = D_del * (D13 * D23 - D12 * D33)
        d_13 = D_del * (D12 * D23 - D13 * D22)
        d_22 = D_del * (D11 * D33 - D13 * D13)
        d_23 = D_del * (D12 * D13 - D11 * D23)
        d_33 = D_del * (D11 * D22 - D12 * D12)
        mp.aC[ix, 1] = B11 * d_11 + B12 * d_12 + B13 * d_13
        mp.aC[ix, 2] = B11 * d_12 + B12 * d_22 + B13 * d_23
        mp.aC[ix, 3] = B11 * d_13 + B12 * d_23 + B13 * d_33
        mp.aC[ix, 4] = B21 * d_11 + B22 * d_12 + B23 * d_13
        mp.aC[ix, 5] = B21 * d_12 + B22 * d_22 + B23 * d_23
        mp.aC[ix, 6] = B21 * d_13 + B22 * d_23 + B23 * d_33
        mp.aC[ix, 7] = B31 * d_11 + B32 * d_12 + B33 * d_13
        mp.aC[ix, 8] = B31 * d_12 + B32 * d_22 + B33 * d_23
        mp.aC[ix, 9] = B31 * d_13 + B32 * d_23 + B33 * d_33
        # update CFL conditions
        nid = attr.nid[ix]
        Ks  = attr.Ks[nid]
        Gs  = attr.Gs[nid]
        cdil = sqrt((Ks + T2(1.333333) * Gs) / mp.ρs[ix]) # 4/3 ≈ 1.333333
        mp.cfl[ix] = min(grid.dx / (cdil + abs(mp.vs[ix, 1])), 
                         grid.dy / (cdil + abs(mp.vs[ix, 2])),
                         grid.dz / (cdil + abs(mp.vs[ix, 3]))) 
    end
end

@kernel inbounds = true function aUpdatestatusvl1_OS!(
    grid   ::    DeviceGrid2D{T1, T2},
    mp     ::DeviceParticle2D{T1, T2},
    ΔT     ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # get values from affine matrix
        dF1, dF2 = mp.aC[ix, 1] * ΔT + T2(1.0), mp.aC[ix, 2] * ΔT
        dF4, dF3 = mp.aC[ix, 4] * ΔT + T2(1.0), mp.aC[ix, 3] * ΔT
        # update deformation gradient incremental matrix
        mp.ΔFs[ix, 1], mp.ΔFs[ix, 2] = dF1, dF2
        mp.ΔFs[ix, 3], mp.ΔFs[ix, 4] = dF3, dF4
        # compute ΔJₚ in the current time step
        ΔJ = dF1 * dF4 - dF2 * dF3 
        # map this value from particle to grid cell
        vol = mp.Ω[ix]
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            NiV = mp.Nij[ix, iy] * vol
            if NiV ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                @KAatomic grid.σm[p2n] += NiV * ΔJ
                @KAatomic grid.Ω[p2n] += NiV
            end
        end
    end
end

@kernel inbounds = true function aUpdatestatusvl1_OS!(
    grid   ::    DeviceGrid3D{T1, T2},
    mp     ::DeviceParticle3D{T1, T2},
    ΔT     ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # get values from affine matrix
        dF1, dF2, dF3 = mp.aC[ix, 1] * ΔT + T2(1.0), mp.aC[ix, 2] * ΔT, mp.aC[ix, 3] * ΔT
        dF5, dF4, dF6 = mp.aC[ix, 5] * ΔT + T2(1.0), mp.aC[ix, 4] * ΔT, mp.aC[ix, 6] * ΔT
        dF9, dF8, dF7 = mp.aC[ix, 9] * ΔT + T2(1.0), mp.aC[ix, 8] * ΔT, mp.aC[ix, 7] * ΔT
        # update deformation gradient incremental matrix
        mp.ΔFs[ix, 1], mp.ΔFs[ix, 2], mp.ΔFs[ix, 3] = dF1, dF2, dF3
        mp.ΔFs[ix, 4], mp.ΔFs[ix, 5], mp.ΔFs[ix, 6] = dF4, dF5, dF6
        mp.ΔFs[ix, 7], mp.ΔFs[ix, 8], mp.ΔFs[ix, 9] = dF7, dF8, dF9
        # compute ΔJₚ in the current time step
        ΔJ = dF1 * dF5 * dF9 + dF2 * dF6 * dF7 + dF3 * dF4 * dF8 - 
             dF7 * dF5 * dF3 - dF8 * dF6 * dF1 - dF9 * dF4 * dF2 
        # map this value from particle to grid cell
        vol = mp.Ω[ix]
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            NiV = mp.Nij[ix, iy] * vol
            if NiV ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                @KAatomic grid.σm[p2n] += NiV * ΔJ
                @KAatomic grid.Ω[p2n] += NiV
            end
        end
    end
end

@kernel inbounds = true function aUpdatestatusvl2_OS!(
    grid   ::    DeviceGrid2D{T1, T2},
    mp     ::DeviceParticle2D{T1, T2},
    ΔT     ::T2
) where {T1, T2}
    ix = @index(Global)
    ΔT_1 = inv(ΔT)
    if ix ≤ mp.np
        co = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Nij = mp.Nij[ix, iy]
            if Nij ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                co += Nij * grid.σm[p2n]
            end
        end
        Jr = mp.ΔFs[ix, 1] * mp.ΔFs[ix, 4] - mp.ΔFs[ix, 2] * mp.ΔFs[ix, 3]
        Jc = (co / Jr) ^ T2(0.5)
        mp.ΔFs[ix, 1] *= Jc; mp.ΔFs[ix, 2] *= Jc; mp.ΔFs[ix, 3] *= Jc; mp.ΔFs[ix, 4] *= Jc
        mp.ΔFs[ix, 1] -= T2(1.0); mp.ΔFs[ix, 4] -= T2(1.0)
        dF1 = mp.ΔFs[ix, 1]; dF2 = mp.ΔFs[ix, 2]; dF3 = mp.ΔFs[ix, 3]; dF4 = mp.ΔFs[ix, 4]
        # strain rate (Second Invariant of Strain Rate Tensor)
        dϵxx = dF1 * ΔT_1
        dϵyy = dF4 * ΔT_1
        dϵxy = T2(0.5) * (dF2 + dF3) * ΔT_1
        mp.ϵv[ix] = sqrt(dϵxx * dϵxx + dϵyy * dϵyy + T2(2.0) * dϵxy * dϵxy)
        # compute strain increment 
        mp.Δϵijs[ix, 1] = dF1
        mp.Δϵijs[ix, 2] = dF4
        mp.Δϵijs[ix, 4] = dF2 + dF3
        # update strain tensor
        mp.ϵijs[ix, 1] += dF1
        mp.ϵijs[ix, 2] += dF4
        mp.ϵijs[ix, 4] += dF2 + dF3
        # deformation gradient matrix
        F1 = mp.F[ix, 1]; F2 = mp.F[ix, 2]; F3 = mp.F[ix, 3]; F4 = mp.F[ix, 4]      
        mp.F[ix, 1] = (dF1 + T2(1.0)) * F1 + dF2 * F3
        mp.F[ix, 2] = (dF1 + T2(1.0)) * F2 + dF2 * F4
        mp.F[ix, 3] = (dF4 + T2(1.0)) * F3 + dF3 * F1
        mp.F[ix, 4] = (dF4 + T2(1.0)) * F4 + dF3 * F2
        # update jacobian value and particle volume
        J = mp.F[ix, 1] * mp.F[ix, 4] - mp.F[ix, 2] * mp.F[ix, 3]
        mp.Ω[ix]  = J * mp.Ω0[ix]
        mp.ρs[ix] = mp.ρs0[ix] / J
    end
end

@kernel inbounds = true function aUpdatestatusvl2_OS!(
    grid   ::    DeviceGrid3D{T1, T2},
    mp     ::DeviceParticle3D{T1, T2},
    ΔT     ::T2
) where {T1, T2}
    ix = @index(Global)
    ΔT_1 = inv(ΔT)
    if ix ≤ mp.np
        co = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Nij = mp.Nij[ix, iy]
            if Nij ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                co += Nij * grid.σm[p2n]
            end
        end
        Jr = mp.ΔFs[ix, 1] * mp.ΔFs[ix, 5] * mp.ΔFs[ix, 9] + 
             mp.ΔFs[ix, 2] * mp.ΔFs[ix, 6] * mp.ΔFs[ix, 7] + 
             mp.ΔFs[ix, 3] * mp.ΔFs[ix, 4] * mp.ΔFs[ix, 8] - 
             mp.ΔFs[ix, 7] * mp.ΔFs[ix, 5] * mp.ΔFs[ix, 3] - 
             mp.ΔFs[ix, 8] * mp.ΔFs[ix, 6] * mp.ΔFs[ix, 1] - 
             mp.ΔFs[ix, 9] * mp.ΔFs[ix, 4] * mp.ΔFs[ix, 2]
        Jc = (co / Jr) ^ T2(0.333333)
        mp.ΔFs[ix, 1] *= Jc; mp.ΔFs[ix, 2] *= Jc; mp.ΔFs[ix, 3] *= Jc
        mp.ΔFs[ix, 4] *= Jc; mp.ΔFs[ix, 5] *= Jc; mp.ΔFs[ix, 6] *= Jc
        mp.ΔFs[ix, 7] *= Jc; mp.ΔFs[ix, 8] *= Jc; mp.ΔFs[ix, 9] *= Jc
        mp.ΔFs[ix, 1] -= T2(1.0); mp.ΔFs[ix, 5] -= T2(1.0); mp.ΔFs[ix, 9] -= T2(1.0)
        dF1 = mp.ΔFs[ix, 1]; dF2 = mp.ΔFs[ix, 2]; dF3 = mp.ΔFs[ix, 3]
        dF4 = mp.ΔFs[ix, 4]; dF5 = mp.ΔFs[ix, 5]; dF6 = mp.ΔFs[ix, 6]
        dF7 = mp.ΔFs[ix, 7]; dF8 = mp.ΔFs[ix, 8]; dF9 = mp.ΔFs[ix, 9]
        # strain rate (Second Invariant of Strain Rate Tensor)
        dϵxx = dF1 * ΔT_1
        dϵyy = dF5 * ΔT_1
        dϵzz = dF9 * ΔT_1
        dϵxy = T2(0.5) * (dF2 + dF4) * ΔT_1
        dϵyz = T2(0.5) * (dF6 + dF8) * ΔT_1
        dϵxz = T2(0.5) * (dF3 + dF7) * ΔT_1
        mp.ϵv[ix] = sqrt(dϵxx * dϵxx + dϵyy * dϵyy + dϵzz * dϵzz + 
            T2(2.0) * (dϵxy * dϵxy + dϵyz * dϵyz + dϵxz * dϵxz))
        # update strain increment
        mp.Δϵijs[ix, 1] = dF1
        mp.Δϵijs[ix, 2] = dF5
        mp.Δϵijs[ix, 3] = dF9
        mp.Δϵijs[ix, 4] = dF2 + dF4
        mp.Δϵijs[ix, 5] = dF6 + dF8
        mp.Δϵijs[ix, 6] = dF3 + dF7
        # update strain tensor
        mp.ϵijs[ix, 1] += dF1
        mp.ϵijs[ix, 2] += dF5
        mp.ϵijs[ix, 3] += dF9
        mp.ϵijs[ix, 4] += dF2 + dF4
        mp.ϵijs[ix, 5] += dF6 + dF8
        mp.ϵijs[ix, 6] += dF3 + dF7
        # deformation gradient matrix
        F1, F2, F3 = mp.F[ix, 1], mp.F[ix, 2], mp.F[ix, 3]
        F4, F5, F6 = mp.F[ix, 4], mp.F[ix, 5], mp.F[ix, 6]
        F7, F8, F9 = mp.F[ix, 7], mp.F[ix, 8], mp.F[ix, 9]
        mp.F[ix, 1] = (dF1 + T2(1.0)) * F1 + dF2 * F4 + dF3 * F7
        mp.F[ix, 2] = (dF1 + T2(1.0)) * F2 + dF2 * F5 + dF3 * F8
        mp.F[ix, 3] = (dF1 + T2(1.0)) * F3 + dF2 * F6 + dF3 * F9
        mp.F[ix, 4] = (dF5 + T2(1.0)) * F4 + dF4 * F1 + dF6 * F7
        mp.F[ix, 5] = (dF5 + T2(1.0)) * F5 + dF4 * F2 + dF6 * F8
        mp.F[ix, 6] = (dF5 + T2(1.0)) * F6 + dF4 * F3 + dF6 * F9
        mp.F[ix, 7] = (dF9 + T2(1.0)) * F7 + dF8 * F4 + dF7 * F1
        mp.F[ix, 8] = (dF9 + T2(1.0)) * F8 + dF8 * F5 + dF7 * F2
        mp.F[ix, 9] = (dF9 + T2(1.0)) * F9 + dF8 * F6 + dF7 * F3
        # update jacobian value and particle volume
        J = mp.F[ix, 1] * mp.F[ix, 5] * mp.F[ix, 9] + 
            mp.F[ix, 2] * mp.F[ix, 6] * mp.F[ix, 7] +
            mp.F[ix, 3] * mp.F[ix, 4] * mp.F[ix, 8] - 
            mp.F[ix, 7] * mp.F[ix, 5] * mp.F[ix, 3] -
            mp.F[ix, 8] * mp.F[ix, 6] * mp.F[ix, 1] - 
            mp.F[ix, 9] * mp.F[ix, 4] * mp.F[ix, 2] 
        mp.Ω[ix]  = J * mp.Ω0[ix]
        mp.ρs[ix] = mp.ρs0[ix] / J
    end
end