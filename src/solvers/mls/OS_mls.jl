#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : OS_mls.jl                                                                  |
|  Description: Implementation of the MLS-MPM                                              |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Functions  : 01. resetgridstatus_MLS_OS! [2D & 3D]                                      |
|               02. resetmpstatus_MLS_OS!   [2D & 3D]                                      |
|               03. aP2G_MLS_OS!            [2D & 3D]                                      |
|               04. solvegrid_MLS_OS!       [2D & 3D]                                      |
|               05. aG2P_MLS_OS!            [2D & 3D]                                      |
+==========================================================================================#

export resetgridstatus_MLS_OS!
export resetmpstatus_MLS_OS!
export aP2G_MLS_OS!
export solvegrid_MLS_OS!
export aG2P_MLS_OS!

@kernel inbounds = true function resetgridstatus_MLS_OS!(
    grid::DeviceGrid2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        grid.σm[ix]     = T2(0.0)
        grid.Ω[ix]      = T2(0.0)
        grid.ms[ix]     = T2(0.0)
        grid.vsT[ix, 1] = T2(0.0)
        grid.vsT[ix, 2] = T2(0.0)
    end
end

@kernel inbounds = true function resetgridstatus_MLS_OS!(
    grid::DeviceGrid3D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        grid.σm[ix]     = T2(0.0)
        grid.Ω[ix]      = T2(0.0)
        grid.ms[ix]     = T2(0.0)
        grid.vsT[ix, 1] = T2(0.0)
        grid.vsT[ix, 2] = T2(0.0)
        grid.vsT[ix, 3] = T2(0.0)
    end
end

@kernel inbounds = true function resetmpstatus_MLS_OS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        gdx_1, gdy_1 = inv(grid.dx), inv(grid.dy)
        # update mass and momentum
        mp.ms[ix]    = mp.Ω[ix]  * mp.ρs[ix]
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1]
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2]
        # base index in the grid
        # note the base index needs a shift of 0.5 (see Taichi)
        mξx, mξy = mp.ξ[ix, 1], mp.ξ[ix, 2]
        bnx = unsafe_trunc(T1, fld(mξx - T2(0.5) * grid.dx - grid.x1, grid.dx))
        bny = unsafe_trunc(T1, fld(mξy - T2(0.5) * grid.dy - grid.y1, grid.dy))
        bid = T1(grid.nny * bnx + grid.nny - bny)
        gξx, gξy = grid.ξ[bid, 1], grid.ξ[bid, 2]
        # p2n index
        @KAunroll for iy in 1:9
            # p2n index
            bac = T1(cld(iy, 3) - 1)
            iyc = T1(iy - bac * 3)
            mp.p2n[ix, iy] = T1(bid + bac * grid.nny - (iyc - 1))
        end
        # compute x value in the basis function N(x)
        x1 = gdx_1 * (mξx - gξx)
        x2 = gdx_1 * (mξx - gξx - grid.dx)
        x3 = gdx_1 * (mξx - gξx - grid.dx - grid.dx)
        y1 = gdy_1 * (mξy - gξy)
        y2 = gdy_1 * (mξy - gξy - grid.dy)
        y3 = gdy_1 * (mξy - gξy - grid.dy - grid.dy)
        Nx1, Nx2, Nx3 = bspline2basis(x1, x2, x3)
        Ny1, Ny2, Ny3 = bspline2basis(y1, y2, y3)
        # assign the value (in order)
        mp.Nij[ix, 1], mp.Nij[ix, 2], mp.Nij[ix, 3] = Nx1 * Ny1, Nx1 * Ny2, Nx1 * Ny3
        mp.Nij[ix, 4], mp.Nij[ix, 5], mp.Nij[ix, 6] = Nx2 * Ny1, Nx2 * Ny2, Nx2 * Ny3
        mp.Nij[ix, 7], mp.Nij[ix, 8], mp.Nij[ix, 9] = Nx3 * Ny1, Nx3 * Ny2, Nx3 * Ny3
    end
end

@kernel inbounds = true function resetmpstatus_MLS_OS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        gdx_1, gdy_1, gdz_1 = inv(grid.dx), inv(grid.dy), inv(grid.dz)
        # update mass and momentum
        mp.ms[ix]    = mp.Ω[ix]  * mp.ρs[ix]
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1]
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2]
        mp.ps[ix, 3] = mp.ms[ix] * mp.vs[ix, 3]
        # base index in the grid
        # note the base index needs a shift of 0.5 (see Taichi)
        mξx, mξy, mξz = mp.ξ[ix, 1], mp.ξ[ix, 2], mp.ξ[ix, 3]
        bnx = unsafe_trunc(T1, fld(mξx - T2(0.5) * grid.dx - grid.x1, grid.dx))
        bny = unsafe_trunc(T1, fld(mξy - T2(0.5) * grid.dy - grid.y1, grid.dy))
        bnz = unsafe_trunc(T1, fld(mξz - T2(0.5) * grid.dz - grid.z1, grid.dz))
        bid = T1(grid.nnx * grid.nny * bnz + grid.nny * bnx + bny + 1)
        gξx, gξy, gξz = grid.ξ[bid, 1], grid.ξ[bid, 2], grid.ξ[bid, 3]
        # p2n index
        @KAunroll for k in 0:2
            for j in 0:2
                for i in 0:2
                    iy = 1 + i + 3*j + 9*k
                    mp.p2n[ix, iy] = bid + i + grid.nny*j + grid.nny*grid.nnx*k
                end
            end
        end
        # compute x value in the basis function N(x)
        x1 = gdx_1 * (mξx - gξx)
        x2 = gdx_1 * (mξx - gξx - grid.dx)
        x3 = gdx_1 * (mξx - gξx - grid.dx - grid.dx)
        y1 = gdy_1 * (mξy - gξy)
        y2 = gdy_1 * (mξy - gξy - grid.dy)
        y3 = gdy_1 * (mξy - gξy - grid.dy - grid.dy)
        z1 = gdz_1 * (mξz - gξz)
        z2 = gdz_1 * (mξz - gξz - grid.dz)
        z3 = gdz_1 * (mξz - gξz - grid.dz - grid.dz)
        Nx1, Nx2, Nx3 = bspline2basis(x1, x2, x3)
        Ny1, Ny2, Ny3 = bspline2basis(y1, y2, y3)
        Nz1, Nz2, Nz3 = bspline2basis(z1, z2, z3)
        # assign the value (in order)
        it = Int32(1)
        @KAunroll for k in 1:3 # z-direction
            for i in 1:3       # x-direction
                for j in 1:3   # y-direction
                    mp.Nij[ix, it] = (Nx1, Nx2, Nx3)[i] * 
                                     (Ny1, Ny2, Ny3)[j] * 
                                     (Nz1, Nz2, Nz3)[k]
                    it += Int32(1)
                end
            end
        end
    end
end

@kernel inbounds = true function aP2G_MLS_OS!(
    grid   ::    DeviceGrid2D{T1, T2},
    mp     ::DeviceParticle2D{T1, T2},
    ΔT     ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        gdx_1, gdy_1 = inv(grid.dx), inv(grid.dy)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                NiM = mp.ms[ix] * Ni
                vtx = mp.Ω[ix] * ΔT * T2(-4.0) * gdx_1 * gdx_1
                vty = mp.Ω[ix] * ΔT * T2(-4.0) * gdy_1 * gdy_1
                gx, gy = grid.ξ[p2n, 1], grid.ξ[p2n, 2]
                mx, my = mp.ξ[ix, 1], mp.ξ[ix, 2]
                dx, dy = gx - mx, gy - my
                # compute nodal mass
                @KAatomic grid.ms[p2n] += NiM
                # compute nodal momentum
                aC1 = mp.aC[ix, 1] * mp.ms[ix] + mp.σij[ix, 1] * vtx
                aC2 = mp.aC[ix, 2] * mp.ms[ix] + mp.σij[ix, 4] * vtx
                aC3 = mp.aC[ix, 3] * mp.ms[ix] + mp.σij[ix, 4] * vty
                aC4 = mp.aC[ix, 4] * mp.ms[ix] + mp.σij[ix, 2] * vty
                # vsT here is momentum
                @KAatomic grid.vsT[p2n, 1] += Ni * (mp.ps[ix, 1] + (aC1 * dx + aC2 * dy))
                @KAatomic grid.vsT[p2n, 2] += Ni * (mp.ps[ix, 2] + (aC3 * dx + aC4 * dy))
            end
        end
    end
end

@kernel inbounds = true function aP2G_MLS_OS!(
    grid   ::    DeviceGrid3D{T1, T2},
    mp     ::DeviceParticle3D{T1, T2},
    ΔT     ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        gdx_1, gdy_1, gdz_1 = inv(grid.dx), inv(grid.dy), inv(grid.dz)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                NiM = mp.ms[ix] * Ni
                vtx = mp.Ω[ix] * ΔT * T2(-4.0) * gdx_1 * gdx_1
                vty = mp.Ω[ix] * ΔT * T2(-4.0) * gdy_1 * gdy_1
                vtz = mp.Ω[ix] * ΔT * T2(-4.0) * gdz_1 * gdz_1
                dx = grid.ξ[p2n, 1] - mp.ξ[ix, 1]
                dy = grid.ξ[p2n, 2] - mp.ξ[ix, 2]
                dz = grid.ξ[p2n, 3] - mp.ξ[ix, 3]
                # compute nodal mass
                @KAatomic grid.ms[p2n] += NiM
                # compute nodal momentum
                aC1 = mp.aC[ix, 1] * mp.ms[ix] + mp.σij[ix, 1] * vtx
                aC2 = mp.aC[ix, 2] * mp.ms[ix] + mp.σij[ix, 4] * vtx
                aC3 = mp.aC[ix, 3] * mp.ms[ix] + mp.σij[ix, 6] * vtx
                aC4 = mp.aC[ix, 4] * mp.ms[ix] + mp.σij[ix, 4] * vty
                aC5 = mp.aC[ix, 5] * mp.ms[ix] + mp.σij[ix, 2] * vty
                aC6 = mp.aC[ix, 6] * mp.ms[ix] + mp.σij[ix, 5] * vty
                aC7 = mp.aC[ix, 7] * mp.ms[ix] + mp.σij[ix, 6] * vtz
                aC8 = mp.aC[ix, 8] * mp.ms[ix] + mp.σij[ix, 5] * vtz
                aC9 = mp.aC[ix, 9] * mp.ms[ix] + mp.σij[ix, 3] * vtz
                # vsT here is momentum
                @KAatomic grid.vsT[p2n, 1] += Ni * (mp.ps[ix, 1] + 
                    (aC1 * dx + aC2 * dy + aC3 * dz))
                @KAatomic grid.vsT[p2n, 2] += Ni * (mp.ps[ix, 2] + 
                    (aC4 * dx + aC5 * dy + aC6 * dz))
                @KAatomic grid.vsT[p2n, 3] += Ni * (mp.ps[ix, 3] + 
                    (aC7 * dx + aC8 * dy + aC9 * dz))
            end
        end
    end
end

@kernel inbounds = true function solvegrid_MLS_OS!(
    grid   ::     DeviceGrid2D{T1, T2},
    bc     ::DeviceVBoundary2D{T1, T2},
    gravity::T2,
    ΔT     ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        # update nodal velocity
        grid.vsT[ix, 1] = grid.vsT[ix, 1] * ms_denom 
        grid.vsT[ix, 2] = grid.vsT[ix, 2] * ms_denom + gravity * ΔT
        # boundary condition
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 2] = bc.vy_s_val[ix] : nothing
    end
end

@kernel inbounds = true function solvegrid_MLS_OS!(
    grid   ::     DeviceGrid3D{T1, T2},
    bc     ::DeviceVBoundary3D{T1, T2},
    gravity::T2,
    ΔT     ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        # update nodal velocity
        grid.vsT[ix, 1] = grid.vsT[ix, 1] * ms_denom 
        grid.vsT[ix, 2] = grid.vsT[ix, 2] * ms_denom
        grid.vsT[ix, 3] = grid.vsT[ix, 3] * ms_denom + gravity * ΔT
        # boundary condition
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vz_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 3] = bc.vz_s_val[ix] : nothing
    end
end

@kernel inbounds = true function aG2P_MLS_OS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
    attr::  DeviceProperty{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        gdx_1, gdy_1 = inv(grid.dx), inv(grid.dy)
        vx = vy = T2(0.0)
        B11 = B12 = B21 = B22 = T2(0.0)
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
            end
        end
        # update particle velocity
        mp.vs[ix, 1] = vx
        mp.vs[ix, 2] = vy
        # update particle position
        mp.ξ[ix, 1] += ΔT * vx
        mp.ξ[ix, 2] += ΔT * vy
        # update affine matrix C
        d_11 = T2(4.0) * gdx_1 * gdx_1
        d_22 = T2(4.0) * gdy_1 * gdy_1
        mp.aC[ix, 1] = B11 * d_11
        mp.aC[ix, 2] = B12 * d_22
        mp.aC[ix, 3] = B21 * d_11
        mp.aC[ix, 4] = B22 * d_22
        # update CFL conditions
        nid = attr.nid[ix]
        Ks  = attr.Ks[nid]
        Gs  = attr.Gs[nid]
        cdil = sqrt((Ks + T2(1.333333) * Gs) / mp.ρs[ix]) # 4/3 ≈ 1.333333
        mp.cfl[ix] = min(grid.dx / (cdil + abs(mp.vs[ix, 1])), 
                         grid.dy / (cdil + abs(mp.vs[ix, 2]))) 
    end
end

@kernel inbounds = true function aG2P_MLS_OS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
    attr::  DeviceProperty{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        gdx_1, gdy_1, gdz_1 = inv(grid.dx), inv(grid.dy), inv(grid.dz)
        vx = vy = vz = T2(0.0)
        B11 = B12 = B13 = B21 = B22 = B23 = B31 = B32 = B33 = T2(0.0)
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
        d_11 = T2(4.0) * gdx_1 * gdx_1
        d_22 = T2(4.0) * gdy_1 * gdy_1
        d_33 = T2(4.0) * gdz_1 * gdz_1
        mp.aC[ix, 1], mp.aC[ix, 2], mp.aC[ix, 3] = B11 * d_11, B12 * d_22, B13 * d_33
        mp.aC[ix, 4], mp.aC[ix, 5], mp.aC[ix, 6] = B21 * d_11, B22 * d_22, B23 * d_33
        mp.aC[ix, 7], mp.aC[ix, 8], mp.aC[ix, 9] = B31 * d_11, B32 * d_22, B33 * d_33
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