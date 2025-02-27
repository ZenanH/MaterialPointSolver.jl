#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : extras.jl                                                                  |
|  Description: some useful functions can be used in the MPM procedure                     |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Functions  : 01. initmpstatus! [2/3D linear/uGIMP/bspline basis]                        |
+==========================================================================================#

export initmpstatus!

"""
    initmpstatus!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2}, Val{:linear})

Description:
---
This function will setup the particle to node and particle to cell index for 2D MPM solver 
    with linear basis functions.
"""
@kernel inbounds=true function initmpstatus!(
    grid::    DeviceGrid2D{T1, T2}, 
    mp  ::DeviceParticle2D{T1, T2},
        ::Val{:linear}
) where {T1, T2}
    ix = @index(Global)
    mp.p2c[ix] = unsafe_trunc(T1,
        cld(mp.ξ[ix, 2] - grid.y1, grid.dy) +
        fld(mp.ξ[ix, 1] - grid.x1, grid.dx) * grid.ncy)
    @KAunroll for iy in Int32(1):Int32(mp.NIC)
        p2n = getP2N_linear(grid, mp.p2c[ix], iy)
        mp.p2n[ix, iy] = p2n
        # compute distance between particle and related nodes
        Δdx = mp.ξ[ix, 1] - grid.ξ[mp.p2n[ix, iy], 1]
        Δdy = mp.ξ[ix, 2] - grid.ξ[mp.p2n[ix, iy], 2]
        # compute basis function
        Nx, dNx = linearbasis(Δdx, grid.dx)
        Ny, dNy = linearbasis(Δdy, grid.dy)
        mp.Nij[ix, iy] =  Nx * Ny # shape function
        mp.∂Nx[ix, iy] = dNx * Ny # x-gradient shape function
        mp.∂Ny[ix, iy] = dNy * Nx # y-gradient shape function
    end
end

"""
    initmpstatus!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2}, Val{:linear})

Description:
---
This function will setup the particle to node and particle to cell index for 3D MPM solver.
    with linear basis functions.
"""
@kernel inbounds=true function initmpstatus!(
    grid::    DeviceGrid3D{T1, T2}, 
    mp  ::DeviceParticle3D{T1, T2},
        ::Val{:linear}
) where {T1, T2}
    ix = @index(Global)
    # compute particle to cell and particle to node index
    mp.p2c[ix] = unsafe_trunc(T1, 
        cld(mp.ξ[ix, 2] - grid.y1, grid.dy) +
        fld(mp.ξ[ix, 3] - grid.z1, grid.dz) * grid.ncy * grid.ncx +
        fld(mp.ξ[ix, 1] - grid.x1, grid.dx) * grid.ncy)
    @KAunroll for iy in Int32(1):Int32(mp.NIC)
        p2n = getP2N_linear(grid, mp.p2c[ix], iy)
        mp.p2n[ix, iy] = p2n
        # compute distance between particle and related nodes
        Δdx = mp.ξ[ix, 1] - grid.ξ[mp.p2n[ix, iy], 1]
        Δdy = mp.ξ[ix, 2] - grid.ξ[mp.p2n[ix, iy], 2]
        Δdz = mp.ξ[ix, 3] - grid.ξ[mp.p2n[ix, iy], 3]
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

"""
    initmpstatus!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2}, Val{:uGIMP})

Description:
---
This function will setup the particle to node and particle to cell index for 2D MPM solver 
    with uGIMP basis functions.
"""
@kernel inbounds=true function initmpstatus!(
    grid::    DeviceGrid2D{T1, T2}, 
    mp  ::DeviceParticle2D{T1, T2},
        ::Val{:uGIMP}
) where {T1, T2}
    ix = @index(Global)
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

"""
    initmpstatus!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2}, Val{:uGIMP})

Description:
---
This function will setup the particle to node and particle to cell index for 3D MPM solver 
    with uGIMP basis functions.
"""
@kernel inbounds=true function initmpstatus!(
    grid::    DeviceGrid3D{T1, T2}, 
    mp  ::DeviceParticle3D{T1, T2},
        ::Val{:uGIMP}
) where {T1, T2}
    ix = @index(Global)
    mpξ1 = mp.ξ[ix, 1]
    mpξ2 = mp.ξ[ix, 2]
    mpξ3 = mp.ξ[ix, 3]
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
        Δdx = mpξ1 - grid.ξ[p2n, 1]
        Δdy = mpξ2 - grid.ξ[p2n, 2]
        Δdz = mpξ3 - grid.ξ[p2n, 3]
        # compute basis function
        if abs(Δdx) < (grid.dx + T2(0.5) * mp.dx) &&
           abs(Δdy) < (grid.dy + T2(0.5) * mp.dy) &&
           abs(Δdz) < (grid.dz + T2(0.5) * mp.dz)
            Nx, dNx = uGIMPbasis(Δdx, grid.dx, mp.dx)
            Ny, dNy = uGIMPbasis(Δdy, grid.dy, mp.dy)
            Nz, dNz = uGIMPbasis(Δdz, grid.dz, mp.dz)
            mp.Nij[ix, viy] =  Nx * Ny * Nz
            mp.∂Nx[ix, viy] = dNx * Ny * Nz # x-gradient basis function
            mp.∂Ny[ix, viy] = dNy * Nx * Nz # y-gradient basis function
            mp.∂Nz[ix, viy] = dNz * Nx * Ny # z-gradient basis function
            mp.p2n[ix, viy] = p2n
            viy += T1(1)
        end
        viy > mp.NIC && break
    end
end

@kernel inbounds = true function initmpstatus!(
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

@kernel inbounds = true function initmpstatus!!(
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
        if size(mp.∂Nx, 1) == mp.np
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
end

@kernel inbounds = true function initmpstatus!(
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
                    if size(mp.∂Nx, 1) == mp.np
                        mp.∂Nx[ix, it] = dxs[i] * Nys[j] * Nzs[k]
                        mp.∂Ny[ix, it] = Nxs[i] * dys[j] * Nzs[k]
                        mp.∂Nz[ix, it] = Nxs[i] * Nys[j] * dzs[k]
                    end
                    it += Int32(1)
                end
            end
        end
    end
end

@kernel inbounds = true function initmpstatus!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
        ::Val{:bspline3}
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
        @KAunroll for iy in Int32(1):Int32(16)
            # p2n index
            p2n = getP2N_uGIMP(grid, mp.p2c[ix], iy)
            # compute distance between particle and related nodes
            rx = (mp.ξ[ix, 1] - grid.ξ[p2n, 1]) / grid.dx
            ry = (mp.ξ[ix, 2] - grid.ξ[p2n, 2]) / grid.dy
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

@kernel inbounds = true function initmpstatus!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
        ::Val{:bspline3}
) where {T1, T2}
    ix = @index(Global)
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
            rx = (mpξ1 - grid.ξ[p2n, 1]) / grid.dx
            ry = (mpξ2 - grid.ξ[p2n, 2]) / grid.dy
            rz = (mpξ3 - grid.ξ[p2n, 3]) / grid.dz
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
