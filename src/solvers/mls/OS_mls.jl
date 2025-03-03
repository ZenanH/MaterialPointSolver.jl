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
        # note the base index needs a shift of 0.5h (see Taichi)
        mξx, mξy = mp.ξ[ix, 1], mp.ξ[ix, 2]
        bnx = unsafe_trunc(T1, fld(mξx - T2(0.5) * grid.dx - grid.x1, grid.dx))
        bny = unsafe_trunc(T1, fld(mξy - T2(0.5) * grid.dy - grid.y1, grid.dy))
        bid = T1(grid.nny * bnx + grid.nny - bny)
        gξx, gξy = grid.ξ[bid, 1], grid.ξ[bid, 2]
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
        it = Int32(1)
        Nxs, Nys = (Nx1, Nx2, Nx3), (Ny1, Ny2, Ny3)
        @KAunroll for i in Int32(1):Int32(3)  # x-direction
            for j in Int32(1):Int32(3)        # y-direction
                mp.p2n[ix, it] = bid - (j - Int32(1)) + grid.nny * (i - Int32(1))
                mp.Nij[ix, it] = Nxs[i] * Nys[j]
                it += Int32(1)
            end
        end
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
        # note the base index needs a shift of 0.5h (see Taichi)
        mξx, mξy, mξz = mp.ξ[ix, 1], mp.ξ[ix, 2], mp.ξ[ix, 3]
        bnx = unsafe_trunc(T1, fld(mξx - T2(0.5) * grid.dx - grid.x1, grid.dx))
        bny = unsafe_trunc(T1, fld(mξy - T2(0.5) * grid.dy - grid.y1, grid.dy))
        bnz = unsafe_trunc(T1, fld(mξz - T2(0.5) * grid.dz - grid.z1, grid.dz))
        bid = T1(grid.nnx * grid.nny * bnz + grid.nny * bnx + bny + 1)
        gξx, gξy, gξz = grid.ξ[bid, 1], grid.ξ[bid, 2], grid.ξ[bid, 3]
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
        Nxs, Nys, Nzs = (Nx1, Nx2, Nx3), (Ny1, Ny2, Ny3), (Nz1, Nz2, Nz3)
        @KAunroll for k in Int32(1):Int32(3) # z-direction
            for i in Int32(1):Int32(3)       # x-direction
                for j in Int32(1):Int32(3)   # y-direction
                    mp.p2n[ix, it] = bid +                 (j - Int32(1)) + 
                                     grid.nny *            (i - Int32(1)) + 
                                     grid.nny * grid.nnx * (k - Int32(1))
                    mp.Nij[ix, it] = Nxs[i] * Nys[j] * Nzs[k]
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
        vtx = mp.Ω[ix] * ΔT * T2(-4.0) * gdx_1 * gdx_1
        vty = mp.Ω[ix] * ΔT * T2(-4.0) * gdy_1 * gdy_1
        # get particle mass and position
        mps, mξx, mξy = mp.ms[ix], mp.ξ[ix, 1], mp.ξ[ix, 2]
        # get particle momentum
        mppx, mppy = mp.ps[ix, 1], mp.ps[ix, 2]
        # get particle stresses
        σxx, σyy, σxy = mp.σij[ix, 1], mp.σij[ix, 2], mp.σij[ix, 4]
        # compute nodal momentum
        aC1 = mp.aC[ix, 1] * mps + σxx * vtx
        aC2 = mp.aC[ix, 2] * mps + σxy * vtx
        aC3 = mp.aC[ix, 3] * mps + σxy * vty
        aC4 = mp.aC[ix, 4] * mps + σyy * vty
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            gi = mp.p2n[ix, iy]
            dx, dy = grid.ξ[gi, 1] - mξx, grid.ξ[gi, 2] - mξy
            # compute nodal mass
            @KAatomic grid.ms[gi] += mps * Ni
            # vsT here is momentum
            @KAatomic grid.vsT[gi, 1] += Ni * (mppx + (aC1 * dx + aC2 * dy))
            @KAatomic grid.vsT[gi, 2] += Ni * (mppy + (aC3 * dx + aC4 * dy))
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
        vtx = mp.Ω[ix] * ΔT * T2(-4.0) * gdx_1 * gdx_1
        vty = mp.Ω[ix] * ΔT * T2(-4.0) * gdy_1 * gdy_1
        vtz = mp.Ω[ix] * ΔT * T2(-4.0) * gdz_1 * gdz_1
        # get particle mass and position
        mps, mξx, mξy, mξz = mp.ms[ix], mp.ξ[ix, 1], mp.ξ[ix, 2], mp.ξ[ix, 3]
        # get particle momentum
        mppx, mppy, mppz = mp.ps[ix, 1], mp.ps[ix, 2], mp.ps[ix, 3]
        # get particle stresses
        σxx, σyy, σzz = mp.σij[ix, 1], mp.σij[ix, 2], mp.σij[ix, 3]
        σxy, σyz, σzx = mp.σij[ix, 4], mp.σij[ix, 5], mp.σij[ix, 6]
        # compute nodal momentum
        aC1 = mp.aC[ix, 1] * mps + σxx * vtx
        aC2 = mp.aC[ix, 2] * mps + σxy * vtx
        aC3 = mp.aC[ix, 3] * mps + σzx * vtx
        aC4 = mp.aC[ix, 4] * mps + σxy * vty
        aC5 = mp.aC[ix, 5] * mps + σyy * vty
        aC6 = mp.aC[ix, 6] * mps + σyz * vty
        aC7 = mp.aC[ix, 7] * mps + σzx * vtz
        aC8 = mp.aC[ix, 8] * mps + σyz * vtz
        aC9 = mp.aC[ix, 9] * mps + σzz * vtz
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            gi = mp.p2n[ix, iy]
            dx = grid.ξ[gi, 1] - mξx
            dy = grid.ξ[gi, 2] - mξy
            dz = grid.ξ[gi, 3] - mξz
            # compute nodal mass
            @KAatomic grid.ms[gi] += mps * Ni
            # vsT here is momentum
            @KAatomic grid.vsT[gi, 1] += Ni * (mppx + (aC1 * dx + aC2 * dy + aC3 * dz))
            @KAatomic grid.vsT[gi, 2] += Ni * (mppy + (aC4 * dx + aC5 * dy + aC6 * dz))
            @KAatomic grid.vsT[gi, 3] += Ni * (mppz + (aC7 * dx + aC8 * dy + aC9 * dz))
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
        # prepare variables
        vx = vy = T2(0.0)
        B11 = B12 = B21 = B22 = T2(0.0)
        gdx_1, gdy_1 = inv(grid.dx), inv(grid.dy)
        mx, my = mp.ξ[ix, 1], mp.ξ[ix, 2]
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni, gi = mp.Nij[ix, iy], mp.p2n[ix, iy]
            gmx, gmy = grid.ξ[gi, 1] - mx, grid.ξ[gi, 2] - my
            vNx, vNy = Ni * grid.vsT[gi, 1], Ni * grid.vsT[gi, 2]
            vx += vNx; vy += vNy
            # compute B matrix
            B11 += vNx * gmx; B12 += vNx * gmy
            B21 += vNy * gmx; B22 += vNy * gmy
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
        # prepare variables
        vx = vy = vz = T2(0.0)
        B11 = B12 = B13 = B21 = B22 = B23 = B31 = B32 = B33 = T2(0.0)
        gdx_1, gdy_1, gdz_1 = inv(grid.dx), inv(grid.dy), inv(grid.dz)
        mx, my, mz = mp.ξ[ix, 1], mp.ξ[ix, 2], mp.ξ[ix, 3]
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni, gi = mp.Nij[ix, iy], mp.p2n[ix, iy]
            gmx, gmy, gmz = grid.ξ[gi, 1] - mx, grid.ξ[gi, 2] - my, grid.ξ[gi, 3] - mz
            vNx, vNy, vNz = Ni * grid.vsT[gi, 1], Ni * grid.vsT[gi, 2], Ni * grid.vsT[gi, 3]
            vx += vNx; vy += vNy; vz += vNz
            # compute B matrix
            B11 += vNx * gmx; B12 += vNx * gmy; B13 += vNx * gmz
            B21 += vNy * gmx; B22 += vNy * gmy; B23 += vNy * gmz
            B31 += vNz * gmx; B32 += vNz * gmy; B33 += vNz * gmz
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