#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : extras.jl                                                                  |
|  Description: some useful functions can be used in the MPM procedure                     |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Functions  : initmpstatus! [2/3D linear/uGIMP/bspline2/bspline3]                        |
+==========================================================================================#

export initmpstatus!

"""
    initmpstatus!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2}, Val{:linear})

Description:
---
This function will setup the particle to node and particle to cell index for 2D MPM solver 
    with linear basis functions.
"""
@kernel inbounds = true function initmpstatus!(
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
        # base index in the grid
        mξx, mξy = mp.ξ[ix, 1], mp.ξ[ix, 2]
        bnx = unsafe_trunc(T1, fld(mξx - grid.x1, grid.dx))
        bny = unsafe_trunc(T1, fld(mξy - grid.y1, grid.dy))
        bid = T1(grid.nny * bnx + grid.nny - bny)
        gξx, gξy = grid.ξ[bid, 1], grid.ξ[bid, 2]
        # compute x value in the basis function N(x)
        x1, y1 = mξx - gξx, mξy - gξy
        x2, y2 = grid.dx - x1, grid.dy - y1
        Nx1, Nx2, ∂x1, ∂x2 = linearbasis(x1, x2, grid.dx)
        Ny1, Ny2, ∂y1, ∂y2 = linearbasis(y1, y2, grid.dy)
        # assign the value (in order)
        it = Int32(1)
        Nxs, Nys = (Nx1, Nx2), (Ny1, Ny2)
        dxs, dys = (∂x1, ∂x2), (∂y1, ∂y2)
        @KAunroll for i in Int32(1):Int32(2)  # x-direction
            for j in Int32(1):Int32(2)  # y-direction
                mp.p2n[ix, it] = bid - (j - Int32(1)) + (i - Int32(1)) * grid.nny
                mp.Nij[ix, it] = Nxs[i] * Nys[j]
                mp.∂Nx[ix, it] = dxs[i] * Nys[j]
                mp.∂Ny[ix, it] = Nxs[i] * dys[j]
                it += Int32(1)
            end
        end
    end
end

"""
    initmpstatus!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2}, Val{:linear})

Description:
---
This function will setup the particle to node and particle to cell index for 3D MPM solver.
    with linear basis functions.
"""
@kernel inbounds = true function initmpstatus!(
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
        # base index in the grid
        mξx, mξy, mξz = mp.ξ[ix, 1], mp.ξ[ix, 2], mp.ξ[ix, 3]
        bnx = unsafe_trunc(T1, fld(mξx - grid.x1, grid.dx))
        bny = unsafe_trunc(T1, fld(mξy - grid.y1, grid.dy))
        bnz = unsafe_trunc(T1, fld(mξz - grid.z1, grid.dz))
        bid = T1(grid.nnx * grid.nny * bnz + grid.nny * bnx + bny + 1)
        gξx, gξy, gξz = grid.ξ[bid, 1], grid.ξ[bid, 2], grid.ξ[bid, 3]
        # compute x value in the basis function N(x)
        x1, y1, z1 = mξx - gξx, mξy - gξy, mξz - gξz
        x2, y2, z2 = grid.dx - x1, grid.dy - y1, grid.dz - z1
        Nx1, Nx2, ∂x1, ∂x2 = linearbasis(x1, x2, grid.dx)
        Ny1, Ny2, ∂y1, ∂y2 = linearbasis(y1, y2, grid.dy)
        Nz1, Nz2, ∂z1, ∂z2 = linearbasis(z1, z2, grid.dz)
        # assign the value (in order)
        it = Int32(1)
        Nxs, Nys, Nzs = (Nx1, Nx2), (Ny1, Ny2), (Nz1, Nz2)
        dxs, dys, dzs = (∂x1, ∂x2), (∂y1, ∂y2),(∂z1, ∂z2)
        @KAunroll for k in Int32(1):Int32(2) # z-direction
            for i in Int32(1):Int32(2)       # x-direction
                for j in Int32(1):Int32(2)   # y-direction
                    mp.p2n[ix, it] = bid +                 (j - Int32(1)) +
                                     grid.nny *            (i - Int32(1)) +
                                     grid.nnx * grid.nny * (k - Int32(1))
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

"""
    initmpstatus!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2}, Val{:uGIMP})

Description:
---
This function will setup the particle to node and particle to cell index for 2D MPM solver 
    with uGIMP basis functions.
"""
@kernel inbounds = true function initmpstatus!(
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
        # base index in the grid
        # note the base index needs a shift of 0.5 (see Taichi)
        mξx, mξy = mp.ξ[ix, 1], mp.ξ[ix, 2]
        bnx = unsafe_trunc(T1, fld(mξx - T2(0.5) * mp.dx - grid.x1, grid.dx))
        bny = unsafe_trunc(T1, fld(mξy - T2(0.5) * mp.dy - grid.y1, grid.dy))
        bid = T1(grid.nny * bnx + grid.nny - bny)
        gξx, gξy = grid.ξ[bid, 1], grid.ξ[bid, 2]
        # compute x value in the basis function N(x)
        x1 = mξx - gξx
        x2 = mξx - gξx - grid.dx
        x3 = mξx - gξx - grid.dx - grid.dx
        y1 = mξy - gξy
        y2 = mξy - gξy - grid.dy
        y3 = mξy - gξy - grid.dy - grid.dy
        Nx1, Nx2, Nx3, ∂x1, ∂x2, ∂x3 = uGIMPbasis(x1, x2, x3, grid.dx, mp.dx)
        Ny1, Ny2, Ny3, ∂y1, ∂y2, ∂y3 = uGIMPbasis(y1, y2, y3, grid.dy, mp.dy)
        # assign the value (in order)
        it = Int32(1)
        Nxs = (Nx1, Nx2, Nx3)
        Nys = (Ny1, Ny2, Ny3)
        dxs = (∂x1, ∂x2, ∂x3)
        dys = (∂y1, ∂y2, ∂y3)
        @KAunroll for i in Int32(1):Int32(3)  # x-direction
            for j in Int32(1):Int32(3)        # y-direction
                mp.p2n[ix, it] = bid - (j - Int32(1)) + grid.nny * (i - Int32(1))
                mp.Nij[ix, it] = Nxs[i] * Nys[j]
                mp.∂Nx[ix, it] = dxs[i] * Nys[j]
                mp.∂Ny[ix, it] = Nxs[i] * dys[j]
                it += Int32(1)
            end
        end
    end
end

"""
    initmpstatus!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2}, Val{:uGIMP})

Description:
---
This function will setup the particle to node and particle to cell index for 3D MPM solver 
    with uGIMP basis functions.
"""
@kernel inbounds = true function initmpstatus!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
        ::Val{:uGIMP}
) where {T1, T2}
    ix = @index(Global)
    T3 = T2
    if ix ≤ mp.np
        mpms = mp.Ω[ix] * mp.ρs[ix]
        # update particle mass and momentum
        mp.ms[ix]    = mpms
        mp.ps[ix, 1] = mpms * mp.vs[ix, 1]
        mp.ps[ix, 2] = mpms * mp.vs[ix, 2]
        mp.ps[ix, 3] = mpms * mp.vs[ix, 3]
        # base index in the grid
        # note the base index needs a shift of 0.5 (see Taichi)
        mξx, mξy, mξz = mp.ξ[ix, 1], mp.ξ[ix, 2], mp.ξ[ix, 3]
        bnx = unsafe_trunc(T1, fld(mξx - T2(0.5) * mp.dx - grid.x1, grid.dx))
        bny = unsafe_trunc(T1, fld(mξy - T2(0.5) * mp.dy - grid.y1, grid.dy))
        bnz = unsafe_trunc(T1, fld(mξz - T2(0.5) * mp.dz - grid.z1, grid.dz))
        bid = T1(grid.nnx * grid.nny * bnz + grid.nny * bnx + bny + 1)
        gξx, gξy, gξz = grid.ξ[bid, 1], grid.ξ[bid, 2], grid.ξ[bid, 3]
        # compute x value in the basis function N(x)
        x1 = mξx - gξx
        x2 = mξx - gξx - grid.dx
        x3 = mξx - gξx - grid.dx - grid.dx
        y1 = mξy - gξy
        y2 = mξy - gξy - grid.dy
        y3 = mξy - gξy - grid.dy - grid.dy
        z1 = mξz - gξz
        z2 = mξz - gξz - grid.dz
        z3 = mξz - gξz - grid.dz - grid.dz
        Nx1, Nx2, Nx3, dx1, dx2, dx3 = uGIMPbasis(x1, x2, x3, grid.dx, mp.dx)
        Ny1, Ny2, Ny3, dy1, dy2, dy3 = uGIMPbasis(y1, y2, y3, grid.dy, mp.dy)
        Nz1, Nz2, Nz3, dz1, dz2, dz3 = uGIMPbasis(z1, z2, z3, grid.dz, mp.dz)
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
                    mp.p2n[ix, it] = bid +                 (j - Int32(1)) + 
                                     grid.nny *            (i - Int32(1)) + 
                                     grid.nny * grid.nnx * (k - Int32(1))
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
        it = Int32(1)
        Nxs = (Nx1, Nx2, Nx3)
        Nys = (Ny1, Ny2, Ny3)
        dxs = (∂x1, ∂x2, ∂x3)
        dys = (∂y1, ∂y2, ∂y3)
        @KAunroll for i in Int32(1):Int32(3)  # x-direction
            for j in Int32(1):Int32(3)        # y-direction
                mp.p2n[ix, it] = bid - (j - Int32(1)) + grid.nny * (i - Int32(1))
                mp.Nij[ix, it] = Nxs[i] * Nys[j]
                if size(mp.∂Nx, 1) == mp.np
                    mp.∂Nx[ix, it] = dxs[i] * Nys[j]
                    mp.∂Ny[ix, it] = Nxs[i] * dys[j]
                end
                it += Int32(1)
            end
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
        @KAunroll for k in Int32(1):Int32(3) # z-direction
            for i in Int32(1):Int32(3)       # x-direction
                for j in Int32(1):Int32(3)   # y-direction
                    mp.p2n[ix, it] = bid +                 (j - Int32(1)) + 
                                     grid.nny *            (i - Int32(1)) + 
                                     grid.nny * grid.nnx * (k - Int32(1))
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
        gdx_1, gdy_1 = inv(grid.dx), inv(grid.dy)
        # update mass and momentum
        mp.ms[ix]    = mp.Ω[ix]  * mp.ρs[ix]
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1]
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2]
        # base index in the grid
        # note the base index needs a shift of 1.0h (see Taichi)
        mξx, mξy = mp.ξ[ix, 1], mp.ξ[ix, 2]
        bnx = unsafe_trunc(T1, fld(mξx - grid.dx - grid.x1, grid.dx))
        bny = unsafe_trunc(T1, fld(mξy - grid.dy - grid.y1, grid.dy))
        bid = T1(grid.nny * bnx + grid.nny - bny)
        gξx, gξy = grid.ξ[bid, 1], grid.ξ[bid, 2]
        # compute x value in the basis function N(x)
        rx1 = (mξx - gξx                              ) * gdx_1
        rx2 = (mξx - gξx - grid.dx                    ) * gdx_1
        rx3 = (mξx - gξx - grid.dx - grid.dx          ) * gdx_1
        rx4 = (mξx - gξx - grid.dx - grid.dx - grid.dx) * gdx_1
        ry1 = (mξy - gξy                              ) * gdy_1
        ry2 = (mξy - gξy - grid.dy                    ) * gdy_1
        ry3 = (mξy - gξy - grid.dy - grid.dy          ) * gdy_1
        ry4 = (mξy - gξy - grid.dy - grid.dy - grid.dy) * gdy_1
        tx1 = get_type(gξx                              , grid.x1, grid.x2, grid.dx)
        tx2 = get_type(gξx + grid.dx                    , grid.x1, grid.x2, grid.dx)
        tx3 = get_type(gξx + grid.dx + grid.dx          , grid.x1, grid.x2, grid.dx)
        tx4 = get_type(gξx + grid.dx + grid.dx + grid.dx, grid.x1, grid.x2, grid.dx)
        ty1 = get_type(gξy                              , grid.y1, grid.y2, grid.dy)
        ty2 = get_type(gξy + grid.dy                    , grid.y1, grid.y2, grid.dy)
        ty3 = get_type(gξy + grid.dy + grid.dy          , grid.y1, grid.y2, grid.dy)
        ty4 = get_type(gξy + grid.dy + grid.dy + grid.dy, grid.y1, grid.y2, grid.dy)
        Nx1, ∂x1 = bspline3basis(rx1, grid.dx, tx1)
        Nx2, ∂x2 = bspline3basis(rx2, grid.dx, tx2)
        Nx3, ∂x3 = bspline3basis(rx3, grid.dx, tx3)
        Nx4, ∂x4 = bspline3basis(rx4, grid.dx, tx4)
        Ny1, ∂y1 = bspline3basis(ry1, grid.dy, ty1)
        Ny2, ∂y2 = bspline3basis(ry2, grid.dy, ty2)
        Ny3, ∂y3 = bspline3basis(ry3, grid.dy, ty3)
        Ny4, ∂y4 = bspline3basis(ry4, grid.dy, ty4)
        # assign the value (in order)
        it = Int32(1)
        Nxs, Nys = (Nx1, Nx2, Nx3, Nx4), (Ny1, Ny2, Ny3, Ny4)
        dxs, dys = (∂x1, ∂x2, ∂x3, ∂x4), (∂y1, ∂y2, ∂y3, ∂y4)
        @KAunroll for i in Int32(1):Int32(4)  # x-direction
            for j in Int32(1):Int32(4)        # y-direction
                mp.p2n[ix, it] = bid - (j - Int32(1)) + grid.nny * (i - Int32(1))
                mp.Nij[ix, it] = Nxs[i] * Nys[j]
                mp.∂Nx[ix, it] = dxs[i] * Nys[j]
                mp.∂Ny[ix, it] = Nxs[i] * dys[j]
                it += Int32(1)
            end
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
        gdx_1, gdy_1, gdz_1 = inv(grid.dx), inv(grid.dy), inv(grid.dz)
        # update mass and momentum
        mp.ms[ix]    = mp.Ω[ix]  * mp.ρs[ix]
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1]
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2]
        mp.ps[ix, 3] = mp.ms[ix] * mp.vs[ix, 3]
        # base index in the grid
        # note the base index needs a shift of 1h (see Taichi)
        mξx, mξy, mξz = mp.ξ[ix, 1], mp.ξ[ix, 2], mp.ξ[ix, 3]
        bnx = unsafe_trunc(T1, fld(mξx - grid.dx - grid.x1, grid.dx))
        bny = unsafe_trunc(T1, fld(mξy - grid.dy - grid.y1, grid.dy))
        bnz = unsafe_trunc(T1, fld(mξz - grid.dz - grid.z1, grid.dz))
        bid = T1(grid.nnx * grid.nny * bnz + grid.nny * bnx + bny + 1)
        gξx, gξy, gξz = grid.ξ[bid, 1], grid.ξ[bid, 2], grid.ξ[bid, 3]
        # compute x value in the basis function N(x)
        rx1 = (mξx - gξx                              ) * gdx_1
        rx2 = (mξx - gξx - grid.dx                    ) * gdx_1
        rx3 = (mξx - gξx - grid.dx - grid.dx          ) * gdx_1
        rx4 = (mξx - gξx - grid.dx - grid.dx - grid.dx) * gdx_1
        ry1 = (mξy - gξy                              ) * gdy_1
        ry2 = (mξy - gξy - grid.dy                    ) * gdy_1
        ry3 = (mξy - gξy - grid.dy - grid.dy          ) * gdy_1
        ry4 = (mξy - gξy - grid.dy - grid.dy - grid.dy) * gdy_1
        rz1 = (mξz - gξz                              ) * gdz_1
        rz2 = (mξz - gξz - grid.dz                    ) * gdz_1
        rz3 = (mξz - gξz - grid.dz - grid.dz          ) * gdz_1
        rz4 = (mξz - gξz - grid.dz - grid.dz - grid.dz) * gdz_1
        tx1 = get_type(gξx                              , grid.x1, grid.x2, grid.dx)
        tx2 = get_type(gξx + grid.dx                    , grid.x1, grid.x2, grid.dx)
        tx3 = get_type(gξx + grid.dx + grid.dx          , grid.x1, grid.x2, grid.dx)
        tx4 = get_type(gξx + grid.dx + grid.dx + grid.dx, grid.x1, grid.x2, grid.dx)
        ty1 = get_type(gξy                              , grid.y1, grid.y2, grid.dy)
        ty2 = get_type(gξy + grid.dy                    , grid.y1, grid.y2, grid.dy)
        ty3 = get_type(gξy + grid.dy + grid.dy          , grid.y1, grid.y2, grid.dy)
        ty4 = get_type(gξy + grid.dy + grid.dy + grid.dy, grid.y1, grid.y2, grid.dy)
        tz1 = get_type(gξz                              , grid.z1, grid.z2, grid.dz)
        tz2 = get_type(gξz + grid.dz                    , grid.z1, grid.z2, grid.dz)
        tz3 = get_type(gξz + grid.dz + grid.dz          , grid.z1, grid.z2, grid.dz)
        tz4 = get_type(gξz + grid.dz + grid.dz + grid.dz, grid.z1, grid.z2, grid.dz)
        Nx1, ∂x1 = bspline3basis(rx1, grid.dx, tx1)
        Nx2, ∂x2 = bspline3basis(rx2, grid.dx, tx2)
        Nx3, ∂x3 = bspline3basis(rx3, grid.dx, tx3)
        Nx4, ∂x4 = bspline3basis(rx4, grid.dx, tx4)
        Ny1, ∂y1 = bspline3basis(ry1, grid.dy, ty1)
        Ny2, ∂y2 = bspline3basis(ry2, grid.dy, ty2)
        Ny3, ∂y3 = bspline3basis(ry3, grid.dy, ty3)
        Ny4, ∂y4 = bspline3basis(ry4, grid.dy, ty4)
        Nz1, ∂z1 = bspline3basis(rz1, grid.dz, tz1)
        Nz2, ∂z2 = bspline3basis(rz2, grid.dz, tz2)
        Nz3, ∂z3 = bspline3basis(rz3, grid.dz, tz3)
        Nz4, ∂z4 = bspline3basis(rz4, grid.dz, tz4)
        # assign the value (in order)
        it = Int32(1)
        Nxs, Nys, Nzs = (Nx1, Nx2, Nx3, Nx4), (Ny1, Ny2, Ny3, Ny4), (Nz1, Nz2, Nz3, Nz4)
        dxs, dys, dzs = (∂x1, ∂x2, ∂x3, ∂x4), (∂y1, ∂y2, ∂y3, ∂y4), (∂z1, ∂z2, ∂z3, ∂z4)
        @KAunroll for k in Int32(1):Int32(4) # z-direction
            for i in Int32(1):Int32(4)       # x-direction
                for j in Int32(1):Int32(4)   # y-direction
                    mp.p2n[ix, it] = bid +                 (j - Int32(1)) + 
                                     grid.nny *            (i - Int32(1)) + 
                                     grid.nny * grid.nnx * (k - Int32(1))
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