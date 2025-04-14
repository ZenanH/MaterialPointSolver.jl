#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : utils_TS.jl                                                                |
|  Description: Basic computing functions for two-phase single-point MPM                   |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Functions  : 01. resetgridstatus_TS! [2D]                                               |
|               02. resetgridstatus_TS! [3D]                                               |
|               03. resetmpstatus_TS!   [2D,   linear basis]                               |
|               04. resetmpstatus_TS!   [3D,   linear basis]                               |
|               05. resetmpstatus_TS!   [2D,    uGIMP basis]                               |
|               06. resetmpstatus_TS!   [3D,    uGIMP basis]                               |
|               07. resetmpstatus_TS!   [2D, bspline2 basis]                               |
|               08. resetmpstatus_TS!   [3D, bspline2 basis]                               |
|               09. resetmpstatus_TS!   [2D, bspline3 basis]                               |
|               10. resetmpstatus_TS!   [3D, bspline3 basis]                               |
|               11. P2G_TS!             [2D]                                               |
|               12. P2G_TS!             [3D]                                               |
|               13. solvegrid_TS!       [2D]                                               |
|               14. solvegrid_TS!       [3D]                                               |
|               15. doublemapping1_TS!  [2D]                                               |
|               16. doublemapping1_TS!  [3D]                                               |
|               17. doublemapping2_TS!  [2D]                                               |
|               18. doublemapping2_TS!  [3D]                                               |
|               19. doublemapping3_TS!  [2D]                                               |
|               20. doublemapping3_TS!  [3D]                                               |
|               21. G2P_TS!             [2D]                                               |
|               22. G2P_TS!             [3D]                                               |
|               23. vollock1_TS!        [2D]                                               |
|               24. vollock1_TS!        [3D]                                               |
|               25. vollock2_TS!        [2D]                                               |
|               26. vollock2_TS!        [3D]                                               |
+==========================================================================================#

export resetgridstatus_TS! 
export resetmpstatus_TS!
export P2G_TS!
export solvegrid_TS! 
export doublemapping1_TS!
export doublemapping2_TS!
export doublemapping3_TS!
export G2P_TS! 
export vollock1_TS!
export vollock2_TS!



@inline function SWCC(S_min, S_max, Pw, p_ref, λ)
    return S_min + (S_max - S_min) * (1.0 + (Pw/p_ref)^inv(1-λ))^(-λ)
end

let 
    S_min = 0.125
    S_max = 1
    p_ref = 3e3
    λ     = 0.7
    Pw    = 2000
    # 生成压力范围，从 1 Pa 到 10^4 = 10000 Pa，采用对数分布
    Pw_values = [10.0^(x) for x in range(2, stop=5, length=100)]
    # 对应的饱和度
    S_values = [SWCC(S_min, S_max, Pw, p_ref, λ) for Pw in Pw_values]
    fig = Figure(size=(600, 400))
    ax = Axis(fig[1, 1], xscale=log10)
    p1 = lines!(ax, hcat(Pw_values, S_values), color=S_values, colormap=:berlin25, 
        linewidth=4, label="SWCC")
    axislegend(ax, merge=true, padding=(10, 6, 0, 0))
    display(fig)
end


"""
    resetgridstatus_TS!(grid::DeviceGrid2D)
    resetgridstatus_TS!(grid::DeviceGrid3D)

Description: [TS: two-phase single-point MPM]
---
Reset grid variables.
"""
@kernel function resetgridstatus_TS!(
    grid::DeviceGrid2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        if ix ≤ grid.nc
            grid.σm[ix] = T2(0.0)
            grid.σw[ix] = T2(0.0)
            grid.Ω[ix]  = T2(0.0)
        end
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

@kernel function resetgridstatus_TS!(
    grid::DeviceGrid3D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        if ix ≤ grid.nc
            grid.σm[ix] = T2(0.0)
            grid.σw[ix] = T2(0.0)
            grid.Ω[ix]  = T2(0.0)
        end
        grid.ms[ix]    = T2(0.0)
        grid.mi[ix]    = T2(0.0)
        grid.mw[ix]    = T2(0.0)
        grid.ps[ix, 1] = T2(0.0)
        grid.ps[ix, 2] = T2(0.0)
        grid.ps[ix, 3] = T2(0.0)
        grid.pw[ix, 1] = T2(0.0)
        grid.pw[ix, 2] = T2(0.0)
        grid.pw[ix, 3] = T2(0.0)
        grid.fs[ix, 1] = T2(0.0)
        grid.fs[ix, 2] = T2(0.0)
        grid.fs[ix, 3] = T2(0.0)
        grid.fw[ix, 1] = T2(0.0)
        grid.fw[ix, 2] = T2(0.0)
        grid.fw[ix, 3] = T2(0.0)
        grid.fd[ix, 1] = T2(0.0)
        grid.fd[ix, 2] = T2(0.0)
        grid.fd[ix, 3] = T2(0.0)
    end
end

"""
    resetmpstatus_TS!(grid::DeviceGrid2D, mp::DeviceParticle2D, basis_type)
    resetmpstatus_TS!(grid::DeviceGrid3D, mp::DeviceParticle3D, basis_type)

Description: [TS: two-phase single-point MPM]
---
1) Update particle mass and momentum.
2) Get topology between particle and grid.
3) Compute basis function N(x) and gradient ∂N(x) by `basis_type`.
4) `basis_type` can be:
    - `:linear`   : linear basis function
    - `:uGIMP`    : uGIMP basis function
    - `:bspline2` : 2nd-order B-spline basis function
    - `:bspline3` : 3rd-order B-spline basis function
Note that `basis_type` is a valued type, so it should be passed by `Val{:linear}` or `Val{:uGIMP}`.
"""
@kernel function resetmpstatus_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
        ::Val{:linear}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # update momentum and mass
        mp.ms[ix] = mp.Ω[ix] * mp.ρs[ix]
        mp.mw[ix] = mp.Ω[ix] * mp.ρw[ix]
        mp.mi[ix] = mp.Ω[ix] * ((T2(1.0) - mp.n[ix]) * mp.ρs[ix] + mp.n[ix] * mp.ρw[ix])
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2] *            mp.n[ix]
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
            for j in Int32(1):Int32(2)        # y-direction
                mp.p2n[ix, it] = bid - (j - Int32(1)) + (i - Int32(1)) * grid.nny 
                mp.Nij[ix, it] = Nxs[i] * Nys[j]
                mp.∂Nx[ix, it] = dxs[i] * Nys[j]
                mp.∂Ny[ix, it] = Nxs[i] * dys[j]
                it += Int32(1)
            end
        end
    end
end

@kernel function resetmpstatus_TS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
        ::Val{:linear}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # update momentum and mass
        mp.ms[ix] = mp.Ω[ix] * mp.ρs[ix]
        mp.mw[ix] = mp.Ω[ix] * mp.ρw[ix]
        mp.mi[ix] = mp.Ω[ix] * ((T2(1.0) - mp.n[ix]) * mp.ρs[ix] + mp.n[ix] * mp.ρw[ix])
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 3] = mp.ms[ix] * mp.vs[ix, 3] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2] *            mp.n[ix]
        mp.pw[ix, 3] = mp.mw[ix] * mp.vw[ix, 3] *            mp.n[ix]
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

@kernel function resetmpstatus_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
        ::Val{:uGIMP}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # update mass and momentum
        mp.ms[ix] = mp.Ω[ix] * mp.ρs[ix]
        mp.mw[ix] = mp.Ω[ix] * mp.ρw[ix]
        mp.mi[ix] = mp.Ω[ix] * ((T2(1.0) - mp.n[ix]) * mp.ρs[ix] + mp.n[ix] * mp.ρw[ix])
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2] *            mp.n[ix]
        # base index in the grid
        # note the base index needs a shift of 0.5lp (see Taichi)
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
        Nxs, Nys = (Nx1, Nx2, Nx3), (Ny1, Ny2, Ny3)
        dxs, dys = (∂x1, ∂x2, ∂x3), (∂y1, ∂y2, ∂y3)
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

@kernel function resetmpstatus_TS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
        ::Val{:uGIMP}
) where {T1, T2}
    ix = @index(Global)   
    if ix ≤ mp.np
        # update momentum and mass
        mp.ms[ix] = mp.Ω[ix] * mp.ρs[ix]
        mp.mw[ix] = mp.Ω[ix] * mp.ρw[ix]
        mp.mi[ix] = mp.Ω[ix] * ((T2(1.0) - mp.n[ix]) * mp.ρs[ix] + mp.n[ix] * mp.ρw[ix])
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 3] = mp.ms[ix] * mp.vs[ix, 3] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2] *            mp.n[ix]
        mp.pw[ix, 3] = mp.mw[ix] * mp.vw[ix, 3] *            mp.n[ix]
        # base index in the grid
        # note the base index needs a shift of 0.5lp (see Taichi)
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
        Nxs, Nys, Nzs = (Nx1, Nx2, Nx3), (Ny1, Ny2, Ny3), (Nz1, Nz2, Nz3)
        dxs, dys, dzs = (dx1, dx2, dx3), (dy1, dy2, dy3), (dz1, dz2, dz3)
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

@kernel function resetmpstatus_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
        ::Val{:bspline2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        gdx_1, gdy_1 = inv(grid.dx), inv(grid.dy)
        # update mass and momentum
        mp.ms[ix] = mp.Ω[ix] * mp.ρs[ix]
        mp.mw[ix] = mp.Ω[ix] * mp.ρw[ix]
        mp.mi[ix] = mp.Ω[ix] * ((T2(1.0) - mp.n[ix]) * mp.ρs[ix] + mp.n[ix] * mp.ρw[ix])
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2] *            mp.n[ix]
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
        Nx1, Nx2, Nx3, ∂x1, ∂x2, ∂x3 = bspline2basis(x1, x2, x3, gdx_1)
        Ny1, Ny2, Ny3, ∂y1, ∂y2, ∂y3 = bspline2basis(y1, y2, y3, gdy_1)
        # assign the value (in order)
        it = Int32(1)
        Nxs, Nys = (Nx1, Nx2, Nx3), (Ny1, Ny2, Ny3)
        dxs, dys = (∂x1, ∂x2, ∂x3), (∂y1, ∂y2, ∂y3)
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

@kernel function resetmpstatus_TS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
        ::Val{:bspline2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        gdx_1, gdy_1, gdz_1 = inv(grid.dx), inv(grid.dy), inv(grid.dz)
        # update momentum and mass
        mp.ms[ix] = mp.Ω[ix] * mp.ρs[ix]
        mp.mw[ix] = mp.Ω[ix] * mp.ρw[ix]
        mp.mi[ix] = mp.Ω[ix] * ((T2(1.0) - mp.n[ix]) * mp.ρs[ix] + mp.n[ix] * mp.ρw[ix])
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 3] = mp.ms[ix] * mp.vs[ix, 3] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2] *            mp.n[ix]
        mp.pw[ix, 3] = mp.mw[ix] * mp.vw[ix, 3] *            mp.n[ix]
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
        Nx1, Nx2, Nx3, dx1, dx2, dx3 = bspline2basis(x1, x2, x3, gdx_1)
        Ny1, Ny2, Ny3, dy1, dy2, dy3 = bspline2basis(y1, y2, y3, gdy_1)
        Nz1, Nz2, Nz3, dz1, dz2, dz3 = bspline2basis(z1, z2, z3, gdz_1)
        # assign the value (in order)
        it = Int32(1)
        Nxs, Nys, Nzs = (Nx1, Nx2, Nx3), (Ny1, Ny2, Ny3), (Nz1, Nz2, Nz3)
        dxs, dys, dzs = (dx1, dx2, dx3), (dy1, dy2, dy3), (dz1, dz2, dz3)
        @KAunroll for k in Int32(1):Int32(3) # z-direction
            for i in Int32(1):Int32(3)       # x-direction
                for j in Int32(1):Int32(3)   # y-direction
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

@kernel function resetmpstatus_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
        ::Val{:bspline3}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        gdx_1, gdy_1 = inv(grid.dx), inv(grid.dy)
        # update mass and momentum
        mp.ms[ix] = mp.Ω[ix] * mp.ρs[ix]
        mp.mw[ix] = mp.Ω[ix] * mp.ρw[ix]
        mp.mi[ix] = mp.Ω[ix] * ((T2(1.0) - mp.n[ix]) * mp.ρs[ix] + mp.n[ix] * mp.ρw[ix])
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2] *            mp.n[ix]
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

@kernel function resetmpstatus_TS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
        ::Val{:bspline3}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        gdx_1, gdy_1, gdz_1 = inv(grid.dx), inv(grid.dy), inv(grid.dz)
        # update mass and momentum
        mp.ms[ix] = mp.Ω[ix] * mp.ρs[ix]
        mp.mw[ix] = mp.Ω[ix] * mp.ρw[ix]
        mp.mi[ix] = mp.Ω[ix] * ((T2(1.0) - mp.n[ix]) * mp.ρs[ix] + mp.n[ix] * mp.ρw[ix])
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 3] = mp.ms[ix] * mp.vs[ix, 3] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2] *            mp.n[ix]
        mp.pw[ix, 3] = mp.mw[ix] * mp.vw[ix, 3] *            mp.n[ix]
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

"""
    P2G_TS!(grid::DeviceGrid2D, mp::DeviceParticle2D, gravity)
    P2G_TS!(grid::DeviceGrid3D, mp::DeviceParticle3D, gravity)

Description: [TS: two-phase single-point MPM]
---
`P2G` procedure for scattering the mass, momentum, and forces from particles to grid.
"""
@kernel function P2G_TS!(
    grid   ::    DeviceGrid2D{T1, T2},
    mp     ::DeviceParticle2D{T1, T2},
    attr   ::  DeviceProperty{T1, T2},
    gravity::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        vol, n, σw = mp.Ω[ix], mp.n[ix], mp.σw[ix]
        mps, mpw, mpi = mp.ms[ix], mp.mw[ix], mp.mi[ix]
        mppsx, mppsy, mppwx, mppwy = mp.ps[ix, 1], mp.ps[ix, 2], mp.pw[ix, 1], mp.pw[ix, 2]
        mpvsx, mpvsy, mpvwx, mpvwy = mp.vs[ix, 1], mp.vs[ix, 2], mp.vw[ix, 1], mp.vw[ix, 2]
        drag = (n * mpw * T2(9.8)) / attr.k[attr.nid[ix]]
        σxx, σyy, σxy = mp.σij[ix, 1], mp.σij[ix, 2], mp.σij[ix, 4]
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                ∂Nx = mp.∂Nx[ix, iy]
                ∂Ny = mp.∂Ny[ix, iy]
                p2n = mp.p2n[ix, iy]
                # compute nodal mass
                @KAatomic grid.ms[p2n] += Ni * mps * (T2(1.0) - n)
                @KAatomic grid.mi[p2n] += Ni * mpw *            n
                @KAatomic grid.mw[p2n] += Ni * mpw
                # compute nodal momentum
                @KAatomic grid.ps[p2n, 1] += Ni * mppsx
                @KAatomic grid.ps[p2n, 2] += Ni * mppsy
                @KAatomic grid.pw[p2n, 1] += Ni * mppwx
                @KAatomic grid.pw[p2n, 2] += Ni * mppwy
                # compute nodal total force
                @KAatomic grid.fw[p2n, 1] += -vol * ∂Nx * σw
                @KAatomic grid.fw[p2n, 2] += -vol * ∂Ny * σw + Ni * mpw * gravity
                @KAatomic grid.fs[p2n, 1] += -vol * (∂Nx * (σxx + σw) + ∂Ny * σxy)
                @KAatomic grid.fs[p2n, 2] += -vol * (∂Ny * (σyy + σw) + ∂Nx * σxy) +
                                               Ni * mpi * gravity
                @KAatomic grid.fd[p2n, 1] += Ni * drag * (mpvwx - mpvsx)
                @KAatomic grid.fd[p2n, 2] += Ni * drag * (mpvwy - mpvsy)
            end
        end
    end
end

@kernel function P2G_TS!(
    grid   ::    DeviceGrid3D{T1, T2},
    mp     ::DeviceParticle3D{T1, T2},
    attr   ::  DeviceProperty{T1, T2},
    gravity::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        vol, n, σw = mp.Ω[ix], mp.n[ix], mp.σw[ix]
        mps, mpw, mpi = mp.ms[ix], mp.mw[ix], mp.mi[ix]
        mppsx, mppsy, mppsz = mp.ps[ix, 1], mp.ps[ix, 2], mp.ps[ix, 3]
        mppwx, mppwy, mppwz = mp.pw[ix, 1], mp.pw[ix, 2], mp.pw[ix, 3]
        mpvsx, mpvsy, mpvsz = mp.vs[ix, 1], mp.vs[ix, 2], mp.vs[ix, 3]
        mpvwx, mpvwy, mpvwz = mp.vw[ix, 1], mp.vw[ix, 2], mp.vw[ix, 3]
        drag = (n * mpw * T2(9.8)) / attr.k[attr.nid[ix]]
        σxx, σyy, σzz = mp.σij[ix, 1], mp.σij[ix, 2], mp.σij[ix, 3]
        σxy, σyz, σzx = mp.σij[ix, 4], mp.σij[ix, 5], mp.σij[ix, 6]
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                ∂Nx = mp.∂Nx[ix, iy]
                ∂Ny = mp.∂Ny[ix, iy]
                ∂Nz = mp.∂Nz[ix, iy]
                p2n = mp.p2n[ix, iy]
                # compute nodal mass
                @KAatomic grid.ms[p2n] += Ni * mps * (T2(1.0) - n)
                @KAatomic grid.mi[p2n] += Ni * mpw *            n
                @KAatomic grid.mw[p2n] += Ni * mpw
                # compute nodal momentum
                @KAatomic grid.ps[p2n, 1] += Ni * mppsx
                @KAatomic grid.ps[p2n, 2] += Ni * mppsy
                @KAatomic grid.ps[p2n, 3] += Ni * mppsz
                @KAatomic grid.pw[p2n, 1] += Ni * mppwx
                @KAatomic grid.pw[p2n, 2] += Ni * mppwy
                @KAatomic grid.pw[p2n, 3] += Ni * mppwz
                # compute nodal total force
                @KAatomic grid.fw[p2n, 1] += -vol * ∂Nx * σw
                @KAatomic grid.fw[p2n, 2] += -vol * ∂Ny * σw
                @KAatomic grid.fw[p2n, 3] += -vol * ∂Nz * σw + Ni * mpw * gravity
                @KAatomic grid.fs[p2n, 1] += -vol * (∂Nx * (σxx + σw) + 
                                                     ∂Ny *  σxy + 
                                                     ∂Nz *  σzx)
                @KAatomic grid.fs[p2n, 2] += -vol * (∂Ny * (σyy + σw) + 
                                                     ∂Nx *  σxy + 
                                                     ∂Nz *  σyz)
                @KAatomic grid.fs[p2n, 3] += -vol * (∂Nz * (σzz + σw) + 
                                                     ∂Nx *  σzx + 
                                                     ∂Ny *  σyz) + Ni * mpi * gravity
                @KAatomic grid.fd[p2n, 1] += Ni * drag * (mpvwx - mpvsx)
                @KAatomic grid.fd[p2n, 2] += Ni * drag * (mpvwy - mpvsy)
                @KAatomic grid.fd[p2n, 3] += Ni * drag * (mpvwz - mpvsz)
            end
        end
    end
end


"""
    solvegrid_TS!(grid::DeviceGrid2D, bc::DeviceVBoundary2D, ΔT, ζs, ζw)
    solvegrid_TS!(grid::DeviceGrid3D, bc::DeviceVBoundary3D, ΔT, ζs, ζw)

Description: [TS: two-phase single-point MPM]
---
Solve equations on the grid and apply boundary conditions.
"""
@kernel function solvegrid_TS!(
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
        grid.vw[ix, 1] = grid.pw[ix, 1] * mi_denom
        grid.vw[ix, 2] = grid.pw[ix, 2] * mi_denom
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

@kernel function solvegrid_TS!(
    grid::     DeviceGrid3D{T1, T2},
    bc  ::DeviceVBoundary3D{T1, T2},
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
        bc.vz_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 3] = bc.vz_s_val[ix] : nothing
        bc.vx_w_idx[ix] ≠ T1(0) ? grid.pw[ix, 1] = bc.vx_w_val[ix] : nothing
        bc.vy_w_idx[ix] ≠ T1(0) ? grid.pw[ix, 2] = bc.vy_w_val[ix] : nothing
        bc.vz_w_idx[ix] ≠ T1(0) ? grid.pw[ix, 3] = bc.vz_w_val[ix] : nothing
        # compute nodal velocity
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        grid.vs[ix, 3] = grid.ps[ix, 3] * ms_denom
        grid.vw[ix, 1] = grid.pw[ix, 1] * mi_denom
        grid.vw[ix, 2] = grid.pw[ix, 2] * mi_denom
        grid.vw[ix, 3] = grid.pw[ix, 3] * mi_denom
        # compute damping force
        dampvw = -ζw * sqrt( grid.fw[ix, 1] * grid.fw[ix, 1]  + 
                             grid.fw[ix, 2] * grid.fw[ix, 2]  +
                             grid.fw[ix, 3] * grid.fw[ix, 3] )
        dampvs = -ζs * sqrt((grid.fs[ix, 1] - grid.fw[ix, 1]) * 
                            (grid.fs[ix, 1] - grid.fw[ix, 1]) +
                            (grid.fs[ix, 2] - grid.fw[ix, 2]) *
                            (grid.fs[ix, 2] - grid.fw[ix, 2]) +
                            (grid.fs[ix, 3] - grid.fw[ix, 3]) *
                            (grid.fs[ix, 3] - grid.fw[ix, 3]))
        # compute node acceleration
        awx = mw_denom * (grid.fw[ix, 1] + dampvw * sign(grid.vw[ix, 1]) - grid.fd[ix, 1])
        awy = mw_denom * (grid.fw[ix, 2] + dampvw * sign(grid.vw[ix, 2]) - grid.fd[ix, 2])
        awz = mw_denom * (grid.fw[ix, 3] + dampvw * sign(grid.vw[ix, 3]) - grid.fd[ix, 3])
        asx = ms_denom * (-grid.mi[ix] * awx + grid.fs[ix, 1] + 
            dampvw * sign(grid.vw[ix, 1]) + dampvs * sign(grid.vs[ix, 1]))
        asy = ms_denom * (-grid.mi[ix] * awy + grid.fs[ix, 2] +
            dampvw * sign(grid.vw[ix, 2]) + dampvs * sign(grid.vs[ix, 2]))
        asz = ms_denom * (-grid.mi[ix] * awz + grid.fs[ix, 3] +
            dampvw * sign(grid.vw[ix, 3]) + dampvs * sign(grid.vs[ix, 3]))
        # update nodal temp velocity
        grid.vsT[ix, 1] = grid.vs[ix, 1] + asx * ΔT
        grid.vsT[ix, 2] = grid.vs[ix, 2] + asy * ΔT
        grid.vsT[ix, 3] = grid.vs[ix, 3] + asz * ΔT
        grid.vwT[ix, 1] = grid.vw[ix, 1] + awx * ΔT
        grid.vwT[ix, 2] = grid.vw[ix, 2] + awy * ΔT
        grid.vwT[ix, 3] = grid.vw[ix, 3] + awz * ΔT
        # apply boundary condition
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vz_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 3] = bc.vz_s_val[ix] : nothing
        bc.vx_w_idx[ix] ≠ T1(0) ? grid.vwT[ix, 1] = bc.vx_w_val[ix] : nothing
        bc.vy_w_idx[ix] ≠ T1(0) ? grid.vwT[ix, 2] = bc.vy_w_val[ix] : nothing
        bc.vz_w_idx[ix] ≠ T1(0) ? grid.vwT[ix, 3] = bc.vz_w_val[ix] : nothing
        # reset grid momentum
        grid.ps[ix, 1] = T2(0.0)
        grid.ps[ix, 2] = T2(0.0)
        grid.ps[ix, 3] = T2(0.0)
        grid.pw[ix, 1] = T2(0.0)
        grid.pw[ix, 2] = T2(0.0)
        grid.pw[ix, 3] = T2(0.0)
    end
end

"""
    doublemapping1_TS!(grid::DeviceGrid2D, mp::DeviceParticle2D, attr::DeviceProperty, ΔT,
        FLIP, PIC)
    doublemapping1_TS!(grid::DeviceGrid3D, mp::DeviceParticle3D, attr::DeviceProperty, ΔT,
        FLIP, PIC)

Description: [TS: two-phase single-point MPM]
---
Mapping results from grid to particles.
"""
@kernel function doublemapping1_TS!(
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
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2] *            mp.n[ix]
        # update CFL conditions
        nid = attr.nid[ix]
        Kw  = attr.Kw[nid]
        Ks  = attr.Ks[nid]
        Gs  = attr.Gs[nid]
        nd  = mp.n[ix] < eps(T2) ? T2(0.0) : inv(mp.n[ix])
        Ec  = Ks + T2(1.333333) * Gs # 4/3 ≈ 1.333333
        Eu  = Ec + Kw * nd
        βs  = sqrt((mp.n[ix] * Ec / Kw) / (T2(1.0) - mp.n[ix] + mp.n[ix] * Ec / Kw))
        c1  = sqrt(Eu / ((T2(1.0) - mp.n[ix]) * mp.ρs[ix] + mp.n[ix] * mp.ρw[ix]))
        c2  = βs * sqrt(Kw / mp.ρw[ix])
        mp.cfl[ix] = min(grid.dx / c1, grid.dx / c2, grid.dy / c1, grid.dy / c2)
    end
end

@kernel function doublemapping1_TS!(
    grid    ::    DeviceGrid3D{T1, T2},
    mp      ::DeviceParticle3D{T1, T2},
    attr    ::  DeviceProperty{T1, T2},
    ΔT      ::T2,
    FLIP    ::T2,
    PIC     ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        ξsx = ξsy = ξsz = vsx = vsy = vsz = T2(0.0)
        ξwx = ξwy = ξwz = vwx = vwy = vwz = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                ξsx += Ni *  grid.vsT[p2n, 1]
                ξsy += Ni *  grid.vsT[p2n, 2]
                ξsz += Ni *  grid.vsT[p2n, 3]
                ξwx += Ni *  grid.vwT[p2n, 1]
                ξwy += Ni *  grid.vwT[p2n, 2]
                ξwz += Ni *  grid.vwT[p2n, 3]
                vsx += Ni * (grid.vsT[p2n, 1] - grid.vs[p2n, 1])
                vsy += Ni * (grid.vsT[p2n, 2] - grid.vs[p2n, 2])
                vsz += Ni * (grid.vsT[p2n, 3] - grid.vs[p2n, 3])
                vwx += Ni * (grid.vwT[p2n, 1] - grid.vw[p2n, 1])
                vwy += Ni * (grid.vwT[p2n, 2] - grid.vw[p2n, 2])
                vwz += Ni * (grid.vwT[p2n, 3] - grid.vw[p2n, 3])
            end
        end
        # update particle ξition
        mp.ξ[ix, 1] += ΔT * ξsx
        mp.ξ[ix, 2] += ΔT * ξsy
        mp.ξ[ix, 3] += ΔT * ξsz
        # update particle velocity
        mp.vs[ix, 1] = FLIP * (mp.vs[ix, 1] + vsx) + PIC * ξsx
        mp.vs[ix, 2] = FLIP * (mp.vs[ix, 2] + vsy) + PIC * ξsy
        mp.vs[ix, 3] = FLIP * (mp.vs[ix, 3] + vsz) + PIC * ξsz
        mp.vw[ix, 1] = FLIP * (mp.vw[ix, 1] + vwx) + PIC * ξwx
        mp.vw[ix, 2] = FLIP * (mp.vw[ix, 2] + vwy) + PIC * ξwy
        mp.vw[ix, 3] = FLIP * (mp.vw[ix, 3] + vwz) + PIC * ξwz
        # update particle momentum
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 3] = mp.ms[ix] * mp.vs[ix, 3] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2] *            mp.n[ix]
        mp.pw[ix, 3] = mp.mw[ix] * mp.vw[ix, 3] *            mp.n[ix]
        # update CFL conditions
        nid = attr.nid[ix]
        Kw  = attr.Kw[nid]
        Ks  = attr.Ks[nid]
        Gs  = attr.Gs[nid]
        nd  = mp.n[ix] < eps(T2) ? T2(0.0) : inv(mp.n[ix])
        Ec  = Ks + T2(1.333333) * Gs # 4/3 ≈ 1.333333
        Eu  = Ec + Kw * nd
        βs  = sqrt((mp.n[ix] * Ec / Kw) / (T2(1.0) - mp.n[ix] + mp.n[ix] * Ec / Kw))
        c1  = sqrt(Eu / ((T2(1.0) - mp.n[ix]) * mp.ρs[ix] + mp.n[ix] * mp.ρw[ix]))
        c2  = βs * sqrt(Kw / mp.ρw[ix])
        mp.cfl[ix] = min(grid.dx / c1, grid.dx / c2, 
                         grid.dy / c1, grid.dy / c2,
                         grid.dz / c1, grid.dz / c2)
    end
end

"""
    doublemapping2_TS!(grid::DeviceGrid2D, mp::DeviceParticle2D)
    doublemapping2_TS!(grid::DeviceGrid3D, mp::DeviceParticle3D)

Description: [TS: two-phase single-point MPM]
---
Scatter momentum from particles to grid.
"""
@kernel function doublemapping2_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # update grid momentum
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                @KAatomic grid.ps[p2n, 1] += mp.ps[ix, 1] * Ni
                @KAatomic grid.ps[p2n, 2] += mp.ps[ix, 2] * Ni
                @KAatomic grid.pw[p2n, 1] += mp.pw[ix, 1] * Ni
                @KAatomic grid.pw[p2n, 2] += mp.pw[ix, 2] * Ni
            end
        end
    end
end

@kernel function doublemapping2_TS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # update grid momentum
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                @KAatomic grid.ps[p2n, 1] += mp.ps[ix, 1] * Ni
                @KAatomic grid.ps[p2n, 2] += mp.ps[ix, 2] * Ni
                @KAatomic grid.ps[p2n, 3] += mp.ps[ix, 3] * Ni
                @KAatomic grid.pw[p2n, 1] += mp.pw[ix, 1] * Ni
                @KAatomic grid.pw[p2n, 2] += mp.pw[ix, 2] * Ni
                @KAatomic grid.pw[p2n, 3] += mp.pw[ix, 3] * Ni
            end
        end
    end
end

"""
    doublemapping3_TS!(grid::DeviceGrid2D, bc::DeviceVBoundary2D, ΔT)
    doublemapping3_TS!(grid::DeviceGrid3D, bc::DeviceVBoundary3D, ΔT)

Description: [TS: two-phase single-point MPM]
---
Solve equations on grid for the updated grid momentum.
"""
@kernel function doublemapping3_TS!(
    grid::     DeviceGrid2D{T1, T2},
    bc  ::DeviceVBoundary2D{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni 
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        mi_denom = grid.mi[ix] < eps(T2) ? T2(0.0) : inv(grid.mi[ix])
        mw_denom = grid.mw[ix] < eps(T2) ? T2(0.0) : inv(grid.mw[ix])
        # compute nodal velocities
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        grid.vw[ix, 1] = grid.pw[ix, 1] * mi_denom
        grid.vw[ix, 2] = grid.pw[ix, 2] * mi_denom
        # fixed Dirichlet nodes
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vx_w_idx[ix] ≠ T1(0) ? grid.vw[ix, 1] = bc.vx_w_val[ix] : nothing
        bc.vy_w_idx[ix] ≠ T1(0) ? grid.vw[ix, 2] = bc.vy_w_val[ix] : nothing
        # compute nodal displacement
        grid.Δus[ix, 1] = grid.vs[ix, 1] * ΔT
        grid.Δus[ix, 2] = grid.vs[ix, 2] * ΔT
        grid.Δuw[ix, 1] = grid.vw[ix, 1] * ΔT
        grid.Δuw[ix, 2] = grid.vw[ix, 2] * ΔT
    end
end

@kernel function doublemapping3_TS!(
    grid::     DeviceGrid3D{T1, T2},
    bc  ::DeviceVBoundary3D{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        mi_denom = grid.mi[ix] < eps(T2) ? T2(0.0) : inv(grid.mi[ix])
        mw_denom = grid.mw[ix] < eps(T2) ? T2(0.0) : inv(grid.mw[ix])
        # compute nodal velocities
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        grid.vs[ix, 3] = grid.ps[ix, 3] * ms_denom
        grid.vw[ix, 1] = grid.pw[ix, 1] * mi_denom
        grid.vw[ix, 2] = grid.pw[ix, 2] * mi_denom
        grid.vw[ix, 3] = grid.pw[ix, 3] * mi_denom
        # fixed Dirichlet nodes
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vz_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 3] = bc.vz_s_val[ix] : nothing
        bc.vx_w_idx[ix] ≠ T1(0) ? grid.vw[ix, 1] = bc.vx_w_val[ix] : nothing
        bc.vy_w_idx[ix] ≠ T1(0) ? grid.vw[ix, 2] = bc.vy_w_val[ix] : nothing
        bc.vz_w_idx[ix] ≠ T1(0) ? grid.vw[ix, 3] = bc.vz_w_val[ix] : nothing
        # compute nodal displacement
        grid.Δus[ix, 1] = grid.vs[ix, 1] * ΔT
        grid.Δus[ix, 2] = grid.vs[ix, 2] * ΔT
        grid.Δus[ix, 3] = grid.vs[ix, 3] * ΔT
        grid.Δuw[ix, 1] = grid.vw[ix, 1] * ΔT
        grid.Δuw[ix, 2] = grid.vw[ix, 2] * ΔT
        grid.Δuw[ix, 3] = grid.vw[ix, 3] * ΔT
    end
end

"""
    G2P_TS!(grid::DeviceGrid2D, mp::DeviceParticle2D, attr::DeviceProperty, ΔT)
    G2P_TS!(grid::DeviceGrid3D, mp::DeviceParticle3D, attr::DeviceProperty, ΔT)

Description: [TS: two-phase single-point MPM]
---
Update particle information.
"""
@kernel function G2P_TS!(
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
        # update pore pressure and n
        mp.σw[ix] += (attr.Kw[attr.nid[ix]] / mp.n[ix]) * (
            (T2(1.0) - mp.n[ix]) * (dfs1 + dfs4) + 
                       mp.n[ix]  * (dfw1 + dfw4))
        mp.n[ix] = clamp(T2(1.0) - (T2(1.0) - mp.n[ix]) / ΔJ, T2(0.0), T2(1.0))
    end
end

@kernel function G2P_TS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
    attr::  DeviceProperty{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        ΔT_1 = inv(ΔT)
        dfs1 = dfs2 = dfs3 = dfs4 = dfs5 = dfs6 = dfs7 = dfs8 = dfs9 = T2(0.0)
        dfw1 = dfw2 = dfw3 = dfw4 = dfw5 = dfw6 = dfw7 = dfw8 = dfw9 = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            if mp.Nij[ix, iy] ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                ∂Nx = mp.∂Nx[ix, iy]; ds1 = grid.Δus[p2n, 1]; dw1 = grid.Δuw[p2n, 1]
                ∂Ny = mp.∂Ny[ix, iy]; ds2 = grid.Δus[p2n, 2]; dw2 = grid.Δuw[p2n, 2]
                ∂Nz = mp.∂Nz[ix, iy]; ds3 = grid.Δus[p2n, 3]; dw3 = grid.Δuw[p2n, 3]
                # compute incremental deformation gradient
                dfs1 += ds1 * ∂Nx; dfs2 += ds1 * ∂Ny; dfs3 += ds1 * ∂Nz
                dfs4 += ds2 * ∂Nx; dfs5 += ds2 * ∂Ny; dfs6 += ds2 * ∂Nz
                dfs7 += ds3 * ∂Nx; dfs8 += ds3 * ∂Ny; dfs9 += ds3 * ∂Nz
                dfw1 += dw1 * ∂Nx; dfw2 += dw1 * ∂Ny; dfw3 += dw1 * ∂Nz
                dfw4 += dw2 * ∂Nx; dfw5 += dw2 * ∂Ny; dfw6 += dw2 * ∂Nz
                dfw7 += dw3 * ∂Nx; dfw8 += dw3 * ∂Ny; dfw9 += dw3 * ∂Nz
            end
        end
        mp.ΔFs[ix, 1] = dfs1; mp.Δfs[ix, 2] = dfs2; mp.Δfs[ix, 3] = dfs3
        mp.ΔFs[ix, 4] = dfs4; mp.Δfs[ix, 5] = dfs5; mp.Δfs[ix, 6] = dfs6
        mp.ΔFs[ix, 7] = dfs7; mp.Δfs[ix, 8] = dfs8; mp.Δfs[ix, 9] = dfs9
        # strain rate (Second Invariant of Strain Rate Tensor)
        ϵvxx = dfs1 * ΔT_1
        ϵvyy = dfs5 * ΔT_1
        ϵvzz = dfs9 * ΔT_1
        ϵvxy = T2(0.5) * (dfs2 + dfs4) * ΔT_1
        ϵvyz = T2(0.5) * (dfs6 + dfs8) * ΔT_1
        ϵvxz = T2(0.5) * (dfs3 + dfs7) * ΔT_1
        mp.ϵv[ix] = sqrt(ϵvxx * ϵvxx + ϵvyy * ϵvyy + ϵvzz * ϵvzz + 
            T2(2.0) * (ϵvxy * ϵvxy + ϵvyz * ϵvyz + ϵvxz * ϵvxz))
        # compute strain increment
        mp.Δϵijs[ix, 1] = dfs1
        mp.Δϵijs[ix, 2] = dfs5
        mp.Δϵijs[ix, 3] = dfs9
        mp.Δϵijs[ix, 4] = dfs2 + dfs4
        mp.Δϵijs[ix, 5] = dfs6 + dfs8
        mp.Δϵijs[ix, 6] = dfs3 + dfs7
        mp.Δϵijw[ix, 1] = dfw1
        mp.Δϵijw[ix, 2] = dfw5
        mp.Δϵijw[ix, 3] = dfw9 
        mp.Δϵijw[ix, 4] = dfw2 + dfw4
        mp.Δϵijw[ix, 5] = dfw6 + dfw8
        mp.Δϵijw[ix, 6] = dfw3 + dfw7
        # update strain tensor
        mp.ϵijs[ix, 1] += dfs1
        mp.ϵijs[ix, 2] += dfs5
        mp.ϵijs[ix, 3] += dfs9
        mp.ϵijs[ix, 4] += dfs2 + dfs4
        mp.ϵijs[ix, 5] += dfs6 + dfs8
        mp.ϵijs[ix, 6] += dfs3 + dfs7
        mp.ϵijw[ix, 1] += dfs1
        mp.ϵijw[ix, 2] += dfs5
        mp.ϵijw[ix, 3] += dfs9
        mp.ϵijw[ix, 4] += dfs2 + dfs4
        mp.ϵijw[ix, 5] += dfs6 + dfs8
        mp.ϵijw[ix, 6] += dfs3 + dfs7
        # deformation gradient matrix
        F1 = mp.F[ix, 1]; F2 = mp.F[ix, 2]; F3 = mp.F[ix, 3]
        F4 = mp.F[ix, 4]; F5 = mp.F[ix, 5]; F6 = mp.F[ix, 6]
        F7 = mp.F[ix, 7]; F8 = mp.F[ix, 8]; F9 = mp.F[ix, 9]        
        mp.F[ix, 1] = (dfs1 + T2(1.0)) * F1 + dfs2 * F4 + dfs3 * F7
        mp.F[ix, 2] = (dfs1 + T2(1.0)) * F2 + dfs2 * F5 + dfs3 * F8
        mp.F[ix, 3] = (dfs1 + T2(1.0)) * F3 + dfs2 * F6 + dfs3 * F9
        mp.F[ix, 4] = (dfs5 + T2(1.0)) * F4 + dfs4 * F1 + dfs6 * F7
        mp.F[ix, 5] = (dfs5 + T2(1.0)) * F5 + dfs4 * F2 + dfs6 * F8
        mp.F[ix, 6] = (dfs5 + T2(1.0)) * F6 + dfs4 * F3 + dfs6 * F9
        mp.F[ix, 7] = (dfs9 + T2(1.0)) * F7 + dfs8 * F4 + dfs7 * F1
        mp.F[ix, 8] = (dfs9 + T2(1.0)) * F8 + dfs8 * F5 + dfs7 * F2
        mp.F[ix, 9] = (dfs9 + T2(1.0)) * F9 + dfs8 * F6 + dfs7 * F3
        # update jacobian value and particle Ωume
        J = mp.F[ix, 1] * mp.F[ix, 5] * mp.F[ix, 9] + 
            mp.F[ix, 2] * mp.F[ix, 6] * mp.F[ix, 7] +
            mp.F[ix, 3] * mp.F[ix, 4] * mp.F[ix, 8] - 
            mp.F[ix, 7] * mp.F[ix, 5] * mp.F[ix, 3] -
            mp.F[ix, 8] * mp.F[ix, 6] * mp.F[ix, 1] - 
            mp.F[ix, 9] * mp.F[ix, 4] * mp.F[ix, 2]
        ΔJ        = J * mp.Ω0[ix] / mp.Ω[ix]
        mp.Ω[ix]  = J * mp.Ω0[ix]
        mp.ρs[ix] = mp.ρs0[ix] / J
        mp.ρw[ix] = mp.ρw0[ix] / J
        # update pore pressure and n
        mp.σw[ix] += (attr.Kw[attr.nid[ix]] / mp.n[ix]) * (
            (T2(1.0) - mp.n[ix]) * (dfs1 + dfs5 + dfs9) + 
                       mp.n[ix]  * (dfw1 + dfw5 + dfw9))
        mp.n[ix] = clamp(T2(1.0) - (T2(1.0) - mp.n[ix]) / ΔJ, T2(0.0), T2(1.0))
    end
end

"""
    vollock1_TS!(grid::DeviceGrid2D, mp::DeviceParticle2D)
    vollock1_TS!(grid::DeviceGrid3D, mp::DeviceParticle3D)

Description: [TS: two-phase single-point MPM]
---
Map the stress and pressure from particles to grid.
"""
@kernel function vollock1_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        p2c = grid.nny * floor(T1, (mp.ξ[ix, 1] - grid.x1) / grid.dx) +
                          ceil(T1, (mp.ξ[ix, 2] - grid.y1) / grid.dy)
        vol = mp.Ω[ix]
        @KAatomic grid.σm[p2c] += vol * mp.σm[ix]
        @KAatomic grid.σw[p2c] += vol * mp.σw[ix]
        @KAatomic grid.Ω[p2c]  += vol
    end
end 

@kernel function vollock1_TS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        p2c = floor(T1, (mp.ξ[ix, 3] - grid.z1) / grid.dz) * grid.nnx * grid.nny + 
              floor(T1, (mp.ξ[ix, 1] - grid.x1) / grid.dx) * grid.nny + 
               ceil(T1, (mp.ξ[ix, 2] - grid.y1) / grid.dy)
        vol = mp.Ω[ix]
        @KAatomic grid.σm[p2c] += vol * mp.σm[ix]
        @KAatomic grid.σw[p2c] += vol * mp.σw[ix]
        @KAatomic grid.Ω[p2c]  += vol
    end
end

"""
    vollock2_TS!(grid::DeviceGrid2D, mp::DeviceParticle2D)
    vollock2_TS!(grid::DeviceGrid3D, mp::DeviceParticle3D)

Description: [TS: two-phase single-point MPM]
---
Update the averaged stress tensor on the particles.
"""
@kernel function vollock2_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        p2c = grid.nny * floor(T1, (mp.ξ[ix, 1] - grid.x1) / grid.dx) +
                          ceil(T1, (mp.ξ[ix, 2] - grid.y1) / grid.dy)
        σm = grid.σm[p2c] / grid.Ω[p2c]
        σw = grid.σw[p2c] / grid.Ω[p2c]
        mp.σij[ix, 1] = mp.sij[ix, 1] + σm
        mp.σij[ix, 2] = mp.sij[ix, 2] + σm
        mp.σij[ix, 3] = mp.sij[ix, 3] + σm
        # update mean stress tensor
        σm = (mp.σij[ix, 1] + mp.σij[ix, 2] + mp.σij[ix, 3]) * T2(0.333333)
        mp.σm[ix] = σm
        mp.σw[ix] = σw
        # update deviatoric stress tensor
        mp.sij[ix, 1] = mp.σij[ix, 1] - σm
        mp.sij[ix, 2] = mp.σij[ix, 2] - σm
        mp.sij[ix, 3] = mp.σij[ix, 3] - σm
    end
end

@kernel function vollock2_TS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        p2c = floor(T1, (mp.ξ[ix, 3] - grid.z1) / grid.dz) * grid.nnx * grid.nny + 
              floor(T1, (mp.ξ[ix, 1] - grid.x1) / grid.dx) * grid.nny + 
               ceil(T1, (mp.ξ[ix, 2] - grid.y1) / grid.dy)
        σm = grid.σm[p2c] / grid.Ω[p2c]
        σw = grid.σw[p2c] / grid.Ω[p2c]
        mp.σij[ix, 1] = mp.sij[ix, 1] + σm
        mp.σij[ix, 2] = mp.sij[ix, 2] + σm
        mp.σij[ix, 3] = mp.sij[ix, 3] + σm
        # update mean stress tensor
        σm = (mp.σij[ix, 1] + mp.σij[ix, 2] + mp.σij[ix, 3]) * T2(0.333333)
        mp.σm[ix] = σm
        mp.σw[ix] = σw
        # update deviatoric stress tensor
        mp.sij[ix, 1] = mp.σij[ix, 1] - σm
        mp.sij[ix, 2] = mp.σij[ix, 2] - σm
        mp.sij[ix, 3] = mp.σij[ix, 3] - σm
    end
end