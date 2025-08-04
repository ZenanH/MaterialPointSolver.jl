#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: Basis (shape) funcs                                                        |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, ISTE, Université de Lausanne                                   |
|  Maintainer : Zenan Huo                                                                  |
+==========================================================================================#

export linearbasis, resetmpstatus!

@inline Base.@propagate_inbounds function linearbasis(x1::T2, x2::T2, invh::T2) where T2
    N1 = T2(1.0) - x1 * invh
    d1 = -invh
    N2 = T2(1.0) - x2 * invh
    d2 = invh
    return N1, N2, d1, d2
end

@kernel function resetmpstatus!(
    grid, mpts, basistype::Linear, dim::Dim2D, ϵ::Precision{T1, T2}
) where {T1, T2}
    ix = @index(Global, Linear)
    if ix ≤ mpts.np
        # base index in the grid
        mξx, mξy = mpts.ξ[1, ix], mpts.ξ[2, ix]
        cx = unsafe_trunc(T1, floor((mξx - grid.x1) * grid.inv_h)) # column
        cy = unsafe_trunc(T1, floor((mξy - grid.y1) * grid.inv_h)) # row (自下而上)
        bid = cx * grid.ny + cy + 1
        gξx = grid.x1 + cx * grid.h
        gξy = grid.y1 + cy * grid.h
        # compute x value in the basis function N(x)
        x1, y1 = mξx - gξx, mξy - gξy
        x2, y2 = grid.h - x1, grid.h - y1
        Nx1, Nx2, ∂x1, ∂x2 = linearbasis(x1, x2, grid.inv_h)
        Ny1, Ny2, ∂y1, ∂y2 = linearbasis(y1, y2, grid.inv_h)
        # assign the value (in order)
        it = 1
        Nxs, Nys = (Nx1, Nx2), (Ny1, Ny2)
        dxs, dys = (∂x1, ∂x2), (∂y1, ∂y2)
        @KAunroll for i in 1:2     # x-direction
            for j in 1:2 # y-direction
                mpts.p2n[it, ix] = bid + (j - 1) + (i - 1) * grid.ny
                mpts.Nij[it, ix] = Nxs[i] * Nys[j]
                mpts.∂Nx[it, ix] = dxs[i] * Nys[j]
                mpts.∂Ny[it, ix] = Nxs[i] * dys[j]
                it += 1
            end
        end
    end
end

@kernel function resetmpstatus!(
    grid, mpts, basistype::Linear, dim::Dim3D, ϵ::Precision{T1, T2}
) where {T1, T2}
    ix = @index(Global, Linear)
    if ix ≤ mpts.np
        # base index in the grid
        mξx, mξy, mξz = mpts.ξ[1, ix], mpts.ξ[2, ix], mpts.ξ[3, ix]
        cx = unsafe_trunc(T1, floor((mξx - grid.x1) * grid.inv_h)) # column
        cy = unsafe_trunc(T1, floor((mξy - grid.y1) * grid.inv_h)) # row (自下而上)
        cz = unsafe_trunc(T1, floor((mξz - grid.z1) * grid.inv_h)) # layer
        bid = cz * (grid.nx * grid.ny) + cx * grid.ny + cy + 1
        gξx = grid.x1 + cx * grid.h
        gξy = grid.y1 + cy * grid.h
        gξz = grid.z1 + cz * grid.h
        # compute x value in the basis function N(x)
        x1, y1, z1 = mξx - gξx, mξy - gξy, mξz - gξz
        x2, y2, z2 = grid.h - x1, grid.h - y1, grid.h - z1
        Nx1, Nx2, ∂x1, ∂x2 = linearbasis(x1, x2, grid.inv_h)
        Ny1, Ny2, ∂y1, ∂y2 = linearbasis(y1, y2, grid.inv_h)
        Nz1, Nz2, ∂z1, ∂z2 = linearbasis(z1, z2, grid.inv_h)
        # assign the value (in order)
        it = 1
        Nxs, Nys, Nzs = (Nx1, Nx2), (Ny1, Ny2), (Nz1, Nz2)
        dxs, dys, dzs = (∂x1, ∂x2), (∂y1, ∂y2), (∂z1, ∂z2)
        @KAunroll for k in 1:2 # z-direction
            for i in 1:2       # x-direction
                for j in 1:2   # y-direction
                    mpts.p2n[it, ix] = bid     +           (j - 1) +
                                       grid.ny           * (i - 1) +
                                       grid.nx * grid.ny * (k - 1)
                    mpts.Nij[it, ix] = Nxs[i] * Nys[j] * Nzs[k]
                    mpts.∂Nx[it, ix] = dxs[i] * Nys[j] * Nzs[k]
                    mpts.∂Ny[it, ix] = Nxs[i] * dys[j] * Nzs[k]
                    mpts.∂Nz[it, ix] = Nxs[i] * Nys[j] * dzs[k]
                    it += 1
                end
            end
        end 
    end
end