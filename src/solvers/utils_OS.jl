#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : utils_OS.jl                                                                |
|  Description: Basic computing functions for one-phase single-point MPM, and these        |
|               functions are mainly used for MUSL update scheme.                          |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Functions  : 01. resetgridstatus_OS!     [2D]                                           |
|               02. resetgridstatus_OS!     [3D]                                           |
|               03. resetmpstatus_OS!       [2D,  linear basis]                            |
|               04. resetmpstatus_OS!       [3D,  linear basis]                            |
|               05. resetmpstatus_OS!       [2D,   uGIMP basis]                            |
|               06. resetmpstatus_OS_CPU!   [2D,   uGIMP basis]                            |
|               07. resetmpstatus_OS!       [3D,   uGIMP basis]                            |
|               08. resetmpstatus_OS_CPU!   [3D,   uGIMP basis]                            |
|               09. resetmpstatus_OS_CPU!   [2D, bspline basis]                            |
|               10. resetmpstatus_OS_CPU!   [3D, bspline basis]                            |
|               11. P2G_OS!                 [2D]                                           |
|               12. P2G_OS!                 [3D]                                           |
|               13. solvegrid_OS!           [2D]                                           |
|               14. solvegrid_OS!           [3D]                                           |
|               15. doublemapping1_OS!      [2D]                                           |
|               16. doublemapping1_OS!      [3D]                                           |
|               17. doublemapping2_OS!      [2D]                                           |
|               18. doublemapping2_OS!      [3D]                                           |
|               19. doublemapping3_OS!      [2D]                                           |
|               20. doublemapping3_OS!      [3D]                                           |
|               21. G2P_OS!                 [2D]                                           |
|               22. G2P_OS!                 [3D]                                           |
|               23. G2Pvl1_OS!              [2D]                                           |
|               24. G2Pvl1_OS!              [3D]                                           |
|               25. G2Pvl2_OS!              [2D]                                           |
|               26. G2Pvl2_OS!              [3D]                                           |
|               27. fastdiv!                [- ]                                           |
+==========================================================================================#

export resetgridstatus_OS!
export resetmpstatus_OS!, resetmpstatus_OS_CUDA!
export P2G_OS! 
export solvegrid_OS!
export doublemapping1_OS!
export doublemapping2_OS!
export doublemapping3_OS!
export G2P_OS! 
export G2Pvl1_OS!
export G2Pvl2_OS!
export fastdiv!

"""
    resetgridstatus_OS!(grid::DeviceGrid2D{T1, T2})

Description:
---
Reset some variables for the grid.
"""
@kernel inbounds = true function resetgridstatus_OS!(
    grid::DeviceGrid2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        grid.σm[ix]    = T2(0.0)
        grid.Ω[ix]     = T2(0.0)
        grid.ms[ix]    = T2(0.0)
        grid.ps[ix, 1] = T2(0.0)
        grid.ps[ix, 2] = T2(0.0)
        grid.fs[ix, 1] = T2(0.0)
        grid.fs[ix, 2] = T2(0.0)
    end
end

"""
    resetgridstatus_OS!(grid::DeviceGrid3D{T1, T2})

Description:
---
Reset some variables for the grid.
"""
@kernel inbounds = true function resetgridstatus_OS!(
    grid::DeviceGrid3D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        grid.σm[ix]    = T2(0.0)
        grid.Ω[ix]     = T2(0.0)
        grid.ms[ix]    = T2(0.0)
        grid.ps[ix, 1] = T2(0.0)
        grid.ps[ix, 2] = T2(0.0)
        grid.ps[ix, 3] = T2(0.0)
        grid.fs[ix, 1] = T2(0.0)
        grid.fs[ix, 2] = T2(0.0)
        grid.fs[ix, 3] = T2(0.0)
    end
end

@kernel inbounds = true function resetmpstatus_OS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
        ::Val{:linear}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # update momentum and mass
        mp.ms[ix]    = mp.Ω[ix]  * mp.ρs[ix]
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1]
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2]
        # compute particle to cell and particle to node index
        mp.p2c[ix] = unsafe_trunc(T1,
            cld(mp.ξ[ix, 2] - grid.y1, grid.dy) +
            fld(mp.ξ[ix, 1] - grid.x1, grid.dx) * grid.ncy)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            p2n = getP2N_linear(grid, mp.p2c[ix], iy)
            mp.p2n[ix, iy] = p2n
            # compute distance between particle and related nodes
            Δdx = mp.ξ[ix, 1] - grid.ξ[p2n, 1]
            Δdy = mp.ξ[ix, 2] - grid.ξ[p2n, 2]
            # compute basis function
            Nx, dNx = linearbasis(Δdx, grid.dx)
            Ny, dNy = linearbasis(Δdy, grid.dy)
            mp.Nij[ix, iy] =  Nx * Ny
            mp.∂Nx[ix, iy] = dNx * Ny # x-gradient shape function
            mp.∂Ny[ix, iy] = dNy * Nx # y-gradient shape function
        end
    end
end

@kernel inbounds = true function resetmpstatus_OS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
        ::Val{:linear}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # update momentum and mass
        mp.ms[ix]    = mp.Ω[ix]  * mp.ρs[ix]
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1]
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2]
        mp.ps[ix, 3] = mp.ms[ix] * mp.vs[ix, 3]
        # compute particle to cell and particle to node index
        mp.p2c[ix] = unsafe_trunc(T1, 
            cld(mp.ξ[ix, 2] - grid.y1, grid.dy) +
            fld(mp.ξ[ix, 3] - grid.z1, grid.dz) * grid.ncy * grid.ncx +
            fld(mp.ξ[ix, 1] - grid.x1, grid.dx) * grid.ncy)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            p2n = getP2N_linear(grid, mp.p2c[ix], iy)
            mp.p2n[ix, iy] = p2n
            # compute distance between particle and related nodes
            Δdx = mp.ξ[ix, 1] - grid.ξ[p2n, 1]
            Δdy = mp.ξ[ix, 2] - grid.ξ[p2n, 2]
            Δdz = mp.ξ[ix, 3] - grid.ξ[p2n, 3]
            # compute basis function
            Nx, dNx = linearbasis(Δdx, grid.dx)
            Ny, dNy = linearbasis(Δdy, grid.dy)
            Nz, dNz = linearbasis(Δdz, grid.dz)
            mp.Nij[ix, iy] =  Nx * Ny * Nz
            mp.∂Nx[ix, iy] = dNx * Ny * Nz # x-gradient shape function
            mp.∂Ny[ix, iy] = dNy * Nx * Nz # y-gradient shape function
            mp.∂Nz[ix, iy] = dNz * Nx * Ny # z-gradient shape function
        end
    end
end

"""
    resetmpstatus_OS_CUDA!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2}, 
        ::Val{:uGIMP})

Description:
---
1. Get topology between particle and grid.
2. Compute the value of basis function (uGIMP).
3. Update particle mass and momentum.
"""
@kernel inbounds = true unsafe_indices = true function resetmpstatus_OS_CUDA!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
        ::Val{:uGIMP}
) where {T1, T2}
    N = @uniform prod(@groupsize())
    gI = @index(Group, Linear)
    i = @index(Local, Linear)
    @synchronize
    ix = (gI - 1) * N + i
    # ix = @index(Global)
    if ix ≤ mp.np
        smem = @localmem T2 Int32(13)
        smem[1]  = grid.dx
        smem[2]  = mp.dx
        smem[3]  = T2(1.0) / (T2(4.0) * grid.dx * mp.dx)
        smem[4]  = T2(1.0) / (grid.dx * mp.dx)
        smem[5]  = T2(1.0) / grid.dx
        smem[6]  = T2(0.5) * mp.dx
        smem[7]  = grid.dy
        smem[8]  = mp.dy
        smem[9]  = T2(1.0) / (T2(4.0) * grid.dy * mp.dy)
        smem[10] = T2(1.0) / (grid.dy * mp.dy)
        smem[11] = T2(1.0) / grid.dy
        smem[12] = T2(0.5) * mp.dy
        # update mass and momentum
        mp.ms[ix]    = mp.Ω[ix]  * mp.ρs[ix]
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1]
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2]
        # p2c index
        mp.p2c[ix] = unsafe_trunc(T1,
            cld(mp.ξ[ix, 2] - grid.y1, grid.dy) +
            fld(mp.ξ[ix, 1] - grid.x1, grid.dx) * grid.ncy)
        # reset Nij, ∂Nx, ∂Ny
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            mp.Nij[ix, iy] = T2(0.0)
            mp.∂Nx[ix, iy] = T2(0.0)
            mp.∂Ny[ix, iy] = T2(0.0)
        end
        viy = T1(1)
        @KAunroll for iy in Int32(1):Int32(16)
            # p2n index
            p2n = getP2N_uGIMP(grid, mp.p2c[ix], iy)
            # compute distance between particle and related nodes
            Δdx = mp.ξ[ix, 1] - grid.ξ[p2n, 1]
            Δdy = mp.ξ[ix, 2] - grid.ξ[p2n, 2]
            # compute basis function
            if abs(Δdx) < (grid.dx + T2(0.5) * mp.dx) &&
               abs(Δdy) < (grid.dy + T2(0.5) * mp.dy)
                Nx, dNx = uGIMPbasisx(Δdx, smem)
                Ny, dNy = uGIMPbasisy(Δdy, smem)
                mp.Nij[ix, viy] =  Nx * Ny
                mp.∂Nx[ix, viy] = dNx * Ny # x-gradient shape function
                mp.∂Ny[ix, viy] = dNy * Nx # y-gradient shape function
                mp.p2n[ix, viy] = p2n
                viy += T1(1)
            end
            viy > mp.NIC && break
        end
    end
end

@kernel inbounds = true function resetmpstatus_OS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
        ::Val{:uGIMP}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # update mass and momentum
        mp.ms[ix]    = mp.Ω[ix]  * mp.ρs[ix]
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1]
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2]
        # p2c index
        mp.p2c[ix] = unsafe_trunc(T1,
            cld(mp.ξ[ix, 2] - grid.y1, grid.dy) +
            fld(mp.ξ[ix, 1] - grid.x1, grid.dx) * grid.ncy)
        # reset Nij, ∂Nx, ∂Ny
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            mp.Nij[ix, iy] = T2(0.0)
            mp.∂Nx[ix, iy] = T2(0.0)
            mp.∂Ny[ix, iy] = T2(0.0)
        end
        viy = T1(1)
        @KAunroll for iy in Int32(1):Int32(16)
            # p2n index
            p2n = getP2N_uGIMP(grid, mp.p2c[ix], iy)
            # compute distance between particle and related nodes
            Δdx = mp.ξ[ix, 1] - grid.ξ[p2n, 1]
            Δdy = mp.ξ[ix, 2] - grid.ξ[p2n, 2]
            # compute basis function
            if abs(Δdx) < (grid.dx + T2(0.5) * mp.dx) &&
               abs(Δdy) < (grid.dy + T2(0.5) * mp.dy)
                Nx, dNx = uGIMPbasis(Δdx, grid.dx, mp.dx)
                Ny, dNy = uGIMPbasis(Δdy, grid.dy, mp.dy)
                mp.Nij[ix, viy] =  Nx * Ny
                mp.∂Nx[ix, viy] = dNx * Ny # x-gradient shape function
                mp.∂Ny[ix, viy] = dNy * Nx # y-gradient shape function
                mp.p2n[ix, viy] = p2n
                viy += T1(1)
            end
            viy > mp.NIC && break
        end
    end
end

"""
    resetmpstatus_OS_CUDA!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2},
        ::Val{:uGIMP})

Description:
---
1. Get topology between particle and grid.
2. Compute the value of basis function (uGIMP).
3. Update particle mass and momentum.
"""
@kernel inbounds = true unsafe_indices = true function resetmpstatus_OS_CUDA!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
        ::Val{:uGIMP}
) where {T1, T2}
    T3 = T2
    N = @uniform prod(@groupsize())
    gI = @index(Group, Linear)
    i = @index(Local, Linear)
    @synchronize
    ix = (gI - 1) * N + i
    # ix = @index(Global)
    smem = @localmem T3 Int32(19)
    smem[1]  = T3(grid.dx)
    smem[2]  = T3(mp.dx)
    smem[3]  = T3(1.0 / (4.0 * grid.dx * mp.dx))
    smem[4]  = T3(1.0 / (grid.dx * mp.dx))
    smem[5]  = T3(1.0 / grid.dx)
    smem[6]  = T3(0.5 * mp.dx)
    smem[7]  = T3(grid.dy)
    smem[8]  = T3(mp.dy)
    smem[9]  = T3(1.0 / (4.0 * grid.dy * mp.dy))
    smem[10] = T3(1.0 / (grid.dy * mp.dy))
    smem[11] = T3(1.0 / grid.dy)
    smem[12] = T3(0.5 * mp.dy)
    smem[13] = T3(grid.dz)
    smem[14] = T3(mp.dz)
    smem[15] = T3(1.0 / (4.0 * grid.dz * mp.dz))
    smem[16] = T3(1.0 / (grid.dz* mp.dz))
    smem[17] = T3(1.0 / grid.dz)
    smem[18] = T3(0.5 * mp.dz)
    if ix ≤ mp.np
        mpξ1 = mp.ξ[ix, 1]
        mpξ2 = mp.ξ[ix, 2]
        mpξ3 = mp.ξ[ix, 3]
        mpms = mp.Ω[ix] * mp.ρs[ix]
        # update particle mass and momentum
        mp.ms[ix]    = mpms
        mp.ps[ix, 1] = mpms * mp.vs[ix, 1]
        mp.ps[ix, 2] = mpms * mp.vs[ix, 2]
        mp.ps[ix, 3] = mpms * mp.vs[ix, 3]
        # p2c index
        mp.p2c[ix] = unsafe_trunc(T1,
            cld(mpξ2 - grid.y1, grid.dy) +
            fld(mpξ3 - grid.z1, grid.dz) * grid.ncy * grid.ncx +
            fld(mpξ1 - grid.x1, grid.dx) * grid.ncy)
        # reset Nij, ∂Nx, ∂Ny, ∂Nz
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            mp.Nij[ix, iy] = T2(0.0)
            mp.∂Nx[ix, iy] = T2(0.0)
            mp.∂Ny[ix, iy] = T2(0.0)
            mp.∂Nz[ix, iy] = T2(0.0)
        end
        viy = T1(1)
        @KAunroll for iy in Int32(1):Int32(64)
            # p2n index
            p2n = getP2N_uGIMP(grid, mp.p2c[ix], iy)
            # compute distance between particle and related nodes
            Δdx = T3(mpξ1 - grid.ξ[p2n, 1])
            Δdy = T3(mpξ2 - grid.ξ[p2n, 2])
            Δdz = T3(mpξ3 - grid.ξ[p2n, 3])
            # compute basis function
            if abs(Δdx) < (T3(grid.dx) + T3(0.5) * T3(mp.dx)) &&
               abs(Δdy) < (T3(grid.dy) + T3(0.5) * T3(mp.dy)) &&
               abs(Δdz) < (T3(grid.dz) + T3(0.5) * T3(mp.dz))
                Nx, dNx = uGIMPbasisx(Δdx, smem)
                Ny, dNy = uGIMPbasisy(Δdy, smem)
                Nz, dNz = uGIMPbasisz(Δdz, smem)
                mp.Nij[ix, viy] = T2( Nx * Ny * Nz)
                mp.∂Nx[ix, viy] = T2(dNx * Ny * Nz) # x-gradient basis function
                mp.∂Ny[ix, viy] = T2(dNy * Nx * Nz) # y-gradient basis function
                mp.∂Nz[ix, viy] = T2(dNz * Nx * Ny) # z-gradient basis function
                mp.p2n[ix, viy] = p2n
                viy += T1(1)
            end
            viy > mp.NIC && break
        end
    end
end

@kernel inbounds = true function resetmpstatus_OS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
        ::Val{:uGIMP}
) where {T1, T2}
    ix = @index(Global)
    T3 = T2
    if ix ≤ mp.np
        mpξ1 = mp.ξ[ix, 1]
        mpξ2 = mp.ξ[ix, 2]
        mpξ3 = mp.ξ[ix, 3]
        mpms = mp.Ω[ix] * mp.ρs[ix]
        # update particle mass and momentum
        mp.ms[ix]    = mpms
        mp.ps[ix, 1] = mpms * mp.vs[ix, 1]
        mp.ps[ix, 2] = mpms * mp.vs[ix, 2]
        mp.ps[ix, 3] = mpms * mp.vs[ix, 3]
        # p2c index
        mp.p2c[ix] = unsafe_trunc(T1,
            cld(mpξ2 - grid.y1, grid.dy) +
            fld(mpξ3 - grid.z1, grid.dz) * grid.ncy * grid.ncx +
            fld(mpξ1 - grid.x1, grid.dx) * grid.ncy)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            mp.Nij[ix, iy] = T2(0.0)
            mp.∂Nx[ix, iy] = T2(0.0)
            mp.∂Ny[ix, iy] = T2(0.0)
            mp.∂Nz[ix, iy] = T2(0.0)
        end
        viy = T1(1)
        @KAunroll for iy in Int32(1):Int32(64)
            # p2n index
            p2n = getP2N_uGIMP(grid, mp.p2c[ix], iy)
            # compute distance betwe en particle and related nodes
            Δdx = T3(mpξ1 - grid.ξ[p2n, 1])
            Δdy = T3(mpξ2 - grid.ξ[p2n, 2])
            Δdz = T3(mpξ3 - grid.ξ[p2n, 3])
            # compute basis function
            if abs(Δdx) < (T3(grid.dx) + T3(0.5) * T3(mp.dx)) &&
               abs(Δdy) < (T3(grid.dy) + T3(0.5) * T3(mp.dy)) &&
               abs(Δdz) < (T3(grid.dz) + T3(0.5) * T3(mp.dz))
                Nx, dNx = uGIMPbasis(Δdx, T3(grid.dx), T3(mp.dx))
                Ny, dNy = uGIMPbasis(Δdy, T3(grid.dy), T3(mp.dy))
                Nz, dNz = uGIMPbasis(Δdz, T3(grid.dz), T3(mp.dz))
                mp.Nij[ix, viy] = T2( Nx * Ny * Nz)
                mp.∂Nx[ix, viy] = T2(dNx * Ny * Nz) # x-gradient basis function
                mp.∂Ny[ix, viy] = T2(dNy * Nx * Nz) # y-gradient basis function
                mp.∂Nz[ix, viy] = T2(dNz * Nx * Ny) # z-gradient basis function
                mp.p2n[ix, viy] = p2n
                viy += T1(1)
            end
            viy > mp.NIC && break
        end
    end
end

@kernel inbounds = true function resetmpstatus_OS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
        ::Val{:bspline2}
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
        Nx1, Nx2, Nx3, ∂x1, ∂x2, ∂x3 = bspline2basis(x1, x2, x3, gdx_1)
        Ny1, Ny2, Ny3, ∂y1, ∂y2, ∂y3 = bspline2basis(y1, y2, y3, gdy_1)
        # assign the value (in order)
        mp.Nij[ix, 1], mp.Nij[ix, 2], mp.Nij[ix, 3] = Nx1 * Ny1, Nx1 * Ny2, Nx1 * Ny3
        mp.Nij[ix, 4], mp.Nij[ix, 5], mp.Nij[ix, 6] = Nx2 * Ny1, Nx2 * Ny2, Nx2 * Ny3
        mp.Nij[ix, 7], mp.Nij[ix, 8], mp.Nij[ix, 9] = Nx3 * Ny1, Nx3 * Ny2, Nx3 * Ny3
        # assign the x-gradient (in order)
        mp.∂Nx[ix, 1], mp.∂Nx[ix, 2], mp.∂Nx[ix, 3] = ∂x1 * Ny1, ∂x1 * Ny2, ∂x1 * Ny3
        mp.∂Nx[ix, 4], mp.∂Nx[ix, 5], mp.∂Nx[ix, 6] = ∂x2 * Ny1, ∂x2 * Ny2, ∂x2 * Ny3
        mp.∂Nx[ix, 7], mp.∂Nx[ix, 8], mp.∂Nx[ix, 9] = ∂x3 * Ny1, ∂x3 * Ny2, ∂x3 * Ny3
        # assign the y-gradient (in order)
        mp.∂Ny[ix, 1], mp.∂Ny[ix, 2], mp.∂Ny[ix, 3] = Nx1 * ∂y1, Nx1 * ∂y2, Nx1 * ∂y3
        mp.∂Ny[ix, 4], mp.∂Ny[ix, 5], mp.∂Ny[ix, 6] = Nx2 * ∂y1, Nx2 * ∂y2, Nx2 * ∂y3
        mp.∂Ny[ix, 7], mp.∂Ny[ix, 8], mp.∂Ny[ix, 9] = Nx3 * ∂y1, Nx3 * ∂y2, Nx3 * ∂y3
    end
end

@kernel inbounds = true function resetmpstatus_OS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
        ::Val{:bspline2}
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
                    mp.p2n[ix, iy] = bid + i + grid.nny * j + grid.nny * grid.nnx * k
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
        Nx1, Nx2, Nx3, dx1, dx2, dx3 = bspline2basis(x1, x2, x3, gdx_1)
        Ny1, Ny2, Ny3, dy1, dy2, dy3 = bspline2basis(y1, y2, y3, gdy_1)
        Nz1, Nz2, Nz3, dz1, dz2, dz3 = bspline2basis(z1, z2, z3, gdz_1)
        # assign the value (in order)
        it = Int32(1)
        Nxs = (Nx1, Nx2, Nx3)
        Nys = (Ny1, Ny2, Ny3)
        Nzs = (Nz1, Nz2, Nz3)
        dxs = (dx1, dx2, dx3)
        dys = (dy1, dy2, dy3)
        dzs = (dz1, dz2, dz3)
        @KAunroll for k in 1:3 # z-direction
            for i in 1:3       # x-direction
                for j in 1:3   # y-direction
                    mp.Nij[ix, it] = Nxs[i] * Nys[j] * Nzs[k]
                    mp.∂Nx[ix, it] = dxs[i] * Nys[j] * Nzs[k]
                    mp.∂Ny[ix, it] = Nxs[i] * dys[j] * Nzs[k]
                    mp.∂Nz[ix, it] = Nxs[i] * Nys[j] * dzs[k]
                    it += Int32(1)
                end
            end
        end
    end
end

@kernel inbounds = true function resetmpstatus_OS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
        ::Val{:bspline3}
) where {T1, T2}
    ix = @index(Global)
    T3 = T2
    if ix ≤ mp.np
        # update mass and momentum
        mp.ms[ix]    = mp.Ω[ix]  * mp.ρs[ix]
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1]
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2]
        # p2c index
        mp.p2c[ix] = unsafe_trunc(T1,
            cld(mp.ξ[ix, 2] - grid.y1, grid.dy) +
            fld(mp.ξ[ix, 1] - grid.x1, grid.dx) * grid.ncy)
        @KAunroll for iy in Int32(1):Int32(16)
            # p2n index
            p2n = getP2N_uGIMP(grid, mp.p2c[ix], iy)
            # compute distance between particle and related nodes
            rx = T3((mp.ξ[ix, 1] - grid.ξ[p2n, 1]) / grid.dx)
            ry = T3((mp.ξ[ix, 2] - grid.ξ[p2n, 2]) / grid.dy)
            # compute basis function
            type_vx = get_type(grid.ξ[p2n, 1], grid.x1, grid.x2, grid.dx)
            Nx, dNx = bspline3basis(rx, grid.dx, type_vx)
            type_vy = get_type(grid.ξ[p2n, 2], grid.y1, grid.y2, grid.dy)
            Ny, dNy = bspline3basis(ry, grid.dy, type_vy)
            mp.Nij[ix, iy] =  Nx * Ny
            mp.∂Nx[ix, iy] = dNx * Ny # x-gradient shape function
            mp.∂Ny[ix, iy] = dNy * Nx # y-gradient shape function
            mp.p2n[ix, iy] = p2n
        end
    end
end

@kernel inbounds = true function resetmpstatus_OS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
        ::Val{:bspline3}
) where {T1, T2}
    ix = @index(Global)
    T3 = T2
    if ix ≤ mp.np
        mpξ1 = mp.ξ[ix, 1]
        mpξ2 = mp.ξ[ix, 2]
        mpξ3 = mp.ξ[ix, 3]
        mpms = mp.Ω[ix] * mp.ρs[ix]
        # update particle mass and momentum
        mp.ms[ix]    = mpms
        mp.ps[ix, 1] = mpms * mp.vs[ix, 1]
        mp.ps[ix, 2] = mpms * mp.vs[ix, 2]
        mp.ps[ix, 3] = mpms * mp.vs[ix, 3]
        # p2c index
        mp.p2c[ix] = unsafe_trunc(T1,
            cld(mpξ2 - grid.y1, grid.dy) +
            fld(mpξ3 - grid.z1, grid.dz) * grid.ncy * grid.ncx +
            fld(mpξ1 - grid.x1, grid.dx) * grid.ncy)
        @KAunroll for iy in Int32(1):Int32(64)
            # p2n index
            p2n = getP2N_uGIMP(grid, mp.p2c[ix], iy)
            # compute distance betwe en particle and related nodes
            rx = T3((mpξ1 - grid.ξ[p2n, 1]) / grid.dx)
            ry = T3((mpξ2 - grid.ξ[p2n, 2]) / grid.dy)
            rz = T3((mpξ3 - grid.ξ[p2n, 3]) / grid.dz)
            # compute basis function
            type_vx = get_type(grid.ξ[p2n, 1], grid.x1, grid.x2, grid.dx)
            Nx, dNx = bspline3basis(rx, grid.dx, type_vx)
            type_vy = get_type(grid.ξ[p2n, 2], grid.y1, grid.y2, grid.dy)
            Ny, dNy = bspline3basis(ry, grid.dy, type_vy)
            type_vz = get_type(grid.ξ[p2n, 3], grid.z1, grid.z2, grid.dz)
            Nz, dNz = bspline3basis(rz, grid.dz, type_vz)
            mp.Nij[ix, iy] = T2( Nx * Ny * Nz)
            mp.∂Nx[ix, iy] = T2(dNx * Ny * Nz) # x-gradient basis function
            mp.∂Ny[ix, iy] = T2(dNy * Nx * Nz) # y-gradient basis function
            mp.∂Nz[ix, iy] = T2(dNz * Nx * Ny) # z-gradient basis function
            mp.p2n[ix, iy] = p2n
        end
    end
end

"""
    P2G_OS!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2}, gravity::T2)

Description:
---
P2G procedure for scattering the mass, momentum, and forces from particles to grid.
"""
@kernel inbounds = true function P2G_OS!(
    grid   ::    DeviceGrid2D{T1, T2},
    mp     ::DeviceParticle2D{T1, T2},
    gravity::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                ∂Nx = mp.∂Nx[ix, iy]
                ∂Ny = mp.∂Ny[ix, iy]
                p2n = mp.p2n[ix, iy]
                vol = mp.Ω[ix]
                NiM = mp.ms[ix] * Ni
                # compute nodal mass
                @KAatomic grid.ms[p2n] += NiM
                # compute nodal momentum
                @KAatomic grid.ps[p2n, 1] += Ni * mp.ps[ix, 1]
                @KAatomic grid.ps[p2n, 2] += Ni * mp.ps[ix, 2]
                # compute nodal total force for solid
                @KAatomic grid.fs[p2n, 1] += -vol * (∂Nx * mp.σij[ix, 1]  + 
                                                     ∂Ny * mp.σij[ix, 4])
                @KAatomic grid.fs[p2n, 2] += -vol * (∂Ny * mp.σij[ix, 2]  + 
                                                     ∂Nx * mp.σij[ix, 4]) + NiM * gravity
            end
        end
    end
end

"""
    P2G_OS!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2}, gravity::T2)

Description:
---
P2G procedure for scattering the mass, momentum, and forces from particles to grid.
"""
@kernel inbounds = true function P2G_OS!(
    grid   ::    DeviceGrid3D{T1, T2},
    mp     ::DeviceParticle3D{T1, T2},
    gravity::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                ∂Nx = mp.∂Nx[ix, iy]
                ∂Ny = mp.∂Ny[ix, iy]
                ∂Nz = mp.∂Nz[ix, iy]
                p2n = mp.p2n[ix, iy]
                vol = mp.Ω[ix]
                NiM = mp.ms[ix] * Ni
                # compute nodal mass
                @KAatomic grid.ms[p2n] += NiM
                # compute nodal momentum
                @KAatomic grid.ps[p2n, 1] += Ni * mp.ps[ix, 1]
                @KAatomic grid.ps[p2n, 2] += Ni * mp.ps[ix, 2]
                @KAatomic grid.ps[p2n, 3] += Ni * mp.ps[ix, 3]
                # compute nodal total force for solid
                @KAatomic grid.fs[p2n, 1] += -vol * (∂Nx * mp.σij[ix, 1]  + 
                                                     ∂Ny * mp.σij[ix, 4]  + 
                                                     ∂Nz * mp.σij[ix, 6])
                @KAatomic grid.fs[p2n, 2] += -vol * (∂Ny * mp.σij[ix, 2]  + 
                                                     ∂Nx * mp.σij[ix, 4]  + 
                                                     ∂Nz * mp.σij[ix, 5])
                @KAatomic grid.fs[p2n, 3] += -vol * (∂Nz * mp.σij[ix, 3]  + 
                                                     ∂Nx * mp.σij[ix, 6]  + 
                                                     ∂Ny * mp.σij[ix, 5]) + NiM * gravity
            end
        end
    end
end

"""
    solvegrid_OS!(grid::DeviceGrid2D{T1, T2}, bc::DeviceVBoundary2D{T1, T2}, ΔT::T2, 
        ζs::T2)

Description:
---
1. Solve equations on grid.
2. Add boundary condition.
3. Update particle velocity based on the acceleration.
"""
@kernel inbounds = true function solvegrid_OS!(
    grid::     DeviceGrid2D{T1, T2},
    bc  ::DeviceVBoundary2D{T1, T2},
    ΔT  ::T2,
    ζs  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        # boundary condition
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 2] = bc.vy_s_val[ix] : nothing
        # compute nodal velocity
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        # damping force for solid
        dampvs = -ζs * sqrt(grid.fs[ix, 1] * grid.fs[ix, 1] + 
                            grid.fs[ix, 2] * grid.fs[ix, 2])
        # compute nodal total force for mixture
        Fs_x = grid.fs[ix, 1] + dampvs * sign(grid.vs[ix, 1])
        Fs_y = grid.fs[ix, 2] + dampvs * sign(grid.vs[ix, 2])
        # update nodal velocity
        grid.vsT[ix, 1] = grid.vs[ix, 1] + Fs_x * ΔT * ms_denom
        grid.vsT[ix, 2] = grid.vs[ix, 2] + Fs_y * ΔT * ms_denom
        # boundary condition
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 2] = bc.vy_s_val[ix] : nothing
        # reset grid momentum
        grid.ps[ix, 1] = T2(0.0)
        grid.ps[ix, 2] = T2(0.0)
    end
end

"""
    solvegrid_OS!(grid::DeviceGrid3D{T1, T2}, bc::DeviceVBoundary3D{T1, T2}, ΔT::T2, 
        ζs::T2)

Description:
---
1. Solve equations on grid.
2. Add boundary condition.
3. Update particle velocity based on the acceleration.
"""
@kernel inbounds = true function solvegrid_OS!(
    grid::     DeviceGrid3D{T1, T2},
    bc  ::DeviceVBoundary3D{T1, T2},
    ΔT  ::T2,
    ζs  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        # boundary condition
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vz_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 3] = bc.vz_s_val[ix] : nothing
        # compute nodal velocity
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        grid.vs[ix, 3] = grid.ps[ix, 3] * ms_denom
        # damping force for solid
        dampvs = -ζs * sqrt(grid.fs[ix, 1] * grid.fs[ix, 1] + 
                            grid.fs[ix, 2] * grid.fs[ix, 2] + 
                            grid.fs[ix, 3] * grid.fs[ix, 3])
        # compute nodal total force for mixture
        Fs_x = grid.fs[ix, 1] + dampvs * sign(grid.vs[ix, 1])
        Fs_y = grid.fs[ix, 2] + dampvs * sign(grid.vs[ix, 2])
        Fs_z = grid.fs[ix, 3] + dampvs * sign(grid.vs[ix, 3])
        # update nodal velocity
        grid.vsT[ix, 1] = grid.vs[ix, 1] + Fs_x * ΔT * ms_denom
        grid.vsT[ix, 2] = grid.vs[ix, 2] + Fs_y * ΔT * ms_denom
        grid.vsT[ix, 3] = grid.vs[ix, 3] + Fs_z * ΔT * ms_denom
        # boundary condition
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vz_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 3] = bc.vz_s_val[ix] : nothing
        # reset grid momentum
        grid.ps[ix, 1] = T2(0.0)
        grid.ps[ix, 2] = T2(0.0)
        grid.ps[ix, 3] = T2(0.0)
    end
end

"""
    doublemapping1_OS!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2},
        attr::DeviceProperty{T1, T2}, ΔT::T2, FLIP::T2, PIC::T2)

Description:
---
Mapping results from grid to particles.
"""
@kernel inbounds = true function doublemapping1_OS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
    attr::  DeviceProperty{T1, T2},
    ΔT  ::T2,
    FLIP::T2,
    PIC ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        ξx = ξy = vx = vy = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                ξx += Ni *  grid.vsT[p2n, 1]
                ξy += Ni *  grid.vsT[p2n, 2]
                vx += Ni * (grid.vsT[p2n, 1] - grid.vs[p2n, 1])
                vy += Ni * (grid.vsT[p2n, 2] - grid.vs[p2n, 2])
            end
        end
        # update particle position
        mp.ξ[ix, 1] += ΔT * ξx
        mp.ξ[ix, 2] += ΔT * ξy
        # update particle velocity
        mp.vs[ix, 1] = FLIP * (mp.vs[ix, 1] + vx) + PIC * ξx
        mp.vs[ix, 2] = FLIP * (mp.vs[ix, 2] + vy) + PIC * ξy
        # update particle momentum
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1]
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2]
        # update CFL conditions
        nid = attr.nid[ix]
        Ks  = attr.Ks[nid]
        Gs  = attr.Gs[nid]
        cdil = sqrt((Ks + T2(1.333333) * Gs) / mp.ρs[ix]) # 4/3 ≈ 1.333333
        mp.cfl[ix] = min(grid.dx / (cdil + abs(mp.vs[ix, 1])), 
                         grid.dy / (cdil + abs(mp.vs[ix, 2]))) 
    end
end

"""
    doublemapping1_OS!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2}, 
        attr::DeviceProperty{T1, T2}, ΔT::T2, FLIP::T2, PIC::T2)

Description:
---
Mapping results from grid to particles.
"""
@kernel inbounds = true function doublemapping1_OS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
    attr::  DeviceProperty{T1, T2},
    ΔT  ::T2,
    FLIP::T2,
    PIC ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        ξx = ξy = ξz = vx = vy = vz = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                ξx += Ni *  grid.vsT[p2n, 1]
                ξy += Ni *  grid.vsT[p2n, 2]
                ξz += Ni *  grid.vsT[p2n, 3]
                vx += Ni * (grid.vsT[p2n, 1] - grid.vs[p2n, 1])
                vy += Ni * (grid.vsT[p2n, 2] - grid.vs[p2n, 2])
                vz += Ni * (grid.vsT[p2n, 3] - grid.vs[p2n, 3])
            end
        end
        # update particle position
        mp.ξ[ix, 1] += ΔT * ξx
        mp.ξ[ix, 2] += ΔT * ξy
        mp.ξ[ix, 3] += ΔT * ξz
        # update particle velocity
        mp.vs[ix, 1] = FLIP * (mp.vs[ix, 1] + vx) + PIC * ξx
        mp.vs[ix, 2] = FLIP * (mp.vs[ix, 2] + vy) + PIC * ξy
        mp.vs[ix, 3] = FLIP * (mp.vs[ix, 3] + vz) + PIC * ξz
        # update particle momentum
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1]
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2]
        mp.ps[ix, 3] = mp.ms[ix] * mp.vs[ix, 3]
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

"""
    doublemapping2_OS!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2})

Description:
---
Scatter momentum from particles to grid.
"""
@kernel inbounds = true function doublemapping2_OS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # update particle position & velocity
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                @KAatomic grid.ps[p2n, 1] += mp.ps[ix, 1] * Ni
                @KAatomic grid.ps[p2n, 2] += mp.ps[ix, 2] * Ni
            end
        end
    end
end

"""
    doublemapping2_OS!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2})

Description:
---
Scatter momentum from particles to grid.
"""
@kernel inbounds = true function doublemapping2_OS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # update particle position & velocity
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                @KAatomic grid.ps[p2n, 1] += mp.ps[ix, 1] * Ni
                @KAatomic grid.ps[p2n, 2] += mp.ps[ix, 2] * Ni
                @KAatomic grid.ps[p2n, 3] += mp.ps[ix, 3] * Ni
            end
        end
    end
end

"""
    doublemapping3_OS!(grid::DeviceGrid2D{T1, T2}, bc::DeviceVBoundary2D{T1, T2}, ΔT::T2)

Description:
---
Solve equations on grid.
"""
@kernel inbounds = true function doublemapping3_OS!(
    grid::     DeviceGrid2D{T1, T2},
    bc  ::DeviceVBoundary2D{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        # compute nodal velocities
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        # fixed Dirichlet nodes
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 2] = bc.vy_s_val[ix] : nothing
        # compute nodal displacement
        grid.Δus[ix, 1] = grid.vs[ix, 1] * ΔT
        grid.Δus[ix, 2] = grid.vs[ix, 2] * ΔT
    end
end

"""
    doublemapping3_OS!(grid::DeviceGrid3D{T1, T2}, bc::KernelVBoundary3D{T1, T2}, ΔT::T2)

Description:
---
Solve equations on grid.
"""
@kernel inbounds = true function doublemapping3_OS!(
    grid::     DeviceGrid3D{T1, T2},
    bc  ::DeviceVBoundary3D{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        # compute nodal velocities
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        grid.vs[ix, 3] = grid.ps[ix, 3] * ms_denom
        # fixed Dirichlet nodes
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vz_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 3] = bc.vz_s_val[ix] : nothing
        # compute nodal displacement
        grid.Δus[ix, 1] = grid.vs[ix, 1] * ΔT
        grid.Δus[ix, 2] = grid.vs[ix, 2] * ΔT
        grid.Δus[ix, 3] = grid.vs[ix, 3] * ΔT
    end
end

"""
    G2P_OS!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2})

Description:
---
Update particle information.
"""
@kernel inbounds = true function G2P_OS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        ΔT_1 = inv(ΔT)
        dF1 = dF2 = dF3 = dF4 = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            if mp.Nij[ix, iy] ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                ∂Nx = mp.∂Nx[ix, iy]
                ∂Ny = mp.∂Ny[ix, iy]
                # compute solid incremental deformation gradient
                dF1 += grid.Δus[p2n, 1] * ∂Nx
                dF2 += grid.Δus[p2n, 1] * ∂Ny
                dF3 += grid.Δus[p2n, 2] * ∂Nx
                dF4 += grid.Δus[p2n, 2] * ∂Ny
            end
        end
        mp.ΔFs[ix, 1] = dF1
        mp.ΔFs[ix, 2] = dF2
        mp.ΔFs[ix, 3] = dF3
        mp.ΔFs[ix, 4] = dF4
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

"""
    G2P_OS!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2})

Description:
---
Update particle information.
"""
@kernel inbounds = true function G2P_OS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        ΔT_1 = T2(1.0) / ΔT
        dF1 = dF2 = dF3 = dF4 = dF5 = dF6 = dF7 = dF8 = dF9 = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            if mp.Nij[ix, iy] ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                ∂Nx = mp.∂Nx[ix, iy]; ds1 = grid.Δus[p2n, 1]
                ∂Ny = mp.∂Ny[ix, iy]; ds2 = grid.Δus[p2n, 2]
                ∂Nz = mp.∂Nz[ix, iy]; ds3 = grid.Δus[p2n, 3]
                # compute solid incremental deformation gradient
                dF1 += ds1 * ∂Nx; dF2 += ds1 * ∂Ny; dF3 += ds1 * ∂Nz
                dF4 += ds2 * ∂Nx; dF5 += ds2 * ∂Ny; dF6 += ds2 * ∂Nz
                dF7 += ds3 * ∂Nx; dF8 += ds3 * ∂Ny; dF9 += ds3 * ∂Nz
            end
        end
        mp.ΔFs[ix, 1] = dF1; mp.ΔFs[ix, 2] = dF2; mp.ΔFs[ix, 3] = dF3
        mp.ΔFs[ix, 4] = dF4; mp.ΔFs[ix, 5] = dF5; mp.ΔFs[ix, 6] = dF6
        mp.ΔFs[ix, 7] = dF7; mp.ΔFs[ix, 8] = dF8; mp.ΔFs[ix, 9] = dF9
        # strain rate (Second Invariant of Strain Rate Tensor)
        dϵxx = dF1 * ΔT_1
        dϵyy = dF5 * ΔT_1
        dϵzz = dF9 * ΔT_1
        dϵxy = T2(0.5) * (dF2 + dF4) * ΔT_1
        dϵyz = T2(0.5) * (dF6 + dF8) * ΔT_1
        dϵxz = T2(0.5) * (dF3 + dF7) * ΔT_1
        mp.ϵv[ix] = sqrt(dϵxx * dϵxx + dϵyy * dϵyy + dϵzz * dϵzz + 
            T2(2.0) * (dϵxy * dϵxy + dϵyz * dϵyz + dϵxz * dϵxz))
        # compute strain increment
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
        F1 = mp.F[ix, 1]; F2 = mp.F[ix, 2]; F3 = mp.F[ix, 3]
        F4 = mp.F[ix, 4]; F5 = mp.F[ix, 5]; F6 = mp.F[ix, 6]
        F7 = mp.F[ix, 7]; F8 = mp.F[ix, 8]; F9 = mp.F[ix, 9]        
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

@kernel inbounds = true function G2Pvl1_OS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        dF1 = dF2 = dF3 = dF4 = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            if mp.Nij[ix, iy] ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                ∂Nx = mp.∂Nx[ix, iy]
                ∂Ny = mp.∂Ny[ix, iy]
                # compute solid incremental deformation gradient
                dF1 += grid.Δus[p2n, 1] * ∂Nx
                dF2 += grid.Δus[p2n, 1] * ∂Ny
                dF3 += grid.Δus[p2n, 2] * ∂Nx
                dF4 += grid.Δus[p2n, 2] * ∂Ny
            end
        end
        dF1 += T2(1.0); dF4 += T2(1.0)
        mp.ΔFs[ix, 1] = dF1; mp.ΔFs[ix, 2] = dF2; mp.ΔFs[ix, 3] = dF3; mp.ΔFs[ix, 4] = dF4
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

@kernel inbounds = true function G2Pvl1_OS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        dF1 = dF2 = dF3 = dF4 = dF5 = dF6 = dF7 = dF8 = dF9 = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            if mp.Nij[ix, iy] ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                ∂Nx = mp.∂Nx[ix, iy]; ds1 = grid.Δus[p2n, 1]
                ∂Ny = mp.∂Ny[ix, iy]; ds2 = grid.Δus[p2n, 2]
                ∂Nz = mp.∂Nz[ix, iy]; ds3 = grid.Δus[p2n, 3]
                # compute solid incremental deformation gradient
                dF1 += ds1 * ∂Nx; dF2 += ds1 * ∂Ny; dF3 += ds1 * ∂Nz
                dF4 += ds2 * ∂Nx; dF5 += ds2 * ∂Ny; dF6 += ds2 * ∂Nz
                dF7 += ds3 * ∂Nx; dF8 += ds3 * ∂Ny; dF9 += ds3 * ∂Nz
            end
        end
        dF1 += T2(1.0); dF5 += T2(1.0); dF9 += T2(1.0)
        mp.ΔFs[ix, 1] = dF1; mp.ΔFs[ix, 2] = dF2; mp.ΔFs[ix, 3] = dF3
        mp.ΔFs[ix, 4] = dF4; mp.ΔFs[ix, 5] = dF5; mp.ΔFs[ix, 6] = dF6
        mp.ΔFs[ix, 7] = dF7; mp.ΔFs[ix, 8] = dF8; mp.ΔFs[ix, 9] = dF9
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

@kernel inbounds = true function G2Pvl2_OS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
    ΔT  ::T2
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
        Jc = sign(co) * abs(co) ^ T2(0.5)
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

@kernel inbounds = true function G2Pvl2_OS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
    ΔT  ::T2
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
        Jc = co ^ T2(0.333333)
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
        # compute strain increment
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
        F1 = mp.F[ix, 1]; F2 = mp.F[ix, 2]; F3 = mp.F[ix, 3]
        F4 = mp.F[ix, 4]; F5 = mp.F[ix, 5]; F6 = mp.F[ix, 6]
        F7 = mp.F[ix, 7]; F8 = mp.F[ix, 8]; F9 = mp.F[ix, 9]        
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

@kernel inbounds = true function fastdiv!(grid::DeviceGrid{T1, T2}) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        Ω = grid.Ω[ix] == T2(0.0) ? T2(0.0) : inv(grid.Ω[ix])
        grid.σm[ix] *= Ω
    end
end