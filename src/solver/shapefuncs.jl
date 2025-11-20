export @N∂N, @Nij

@inline Base.@propagate_inbounds function bspline2basis(
    ξ1::T2,
    ξ2::T2, 
    ξ3::T2,
    Δ1::T2
) where T2
    N1 = T2(0.50) * (T2(1.5) - ξ1) * (T2(1.5) - ξ1)
    N2 = T2(0.75) - ξ2 * ξ2
    N3 = T2(0.50) * (T2(1.5) + ξ3) * (T2(1.5) + ξ3)
    d1 = ξ1 - T2(1.5)
    d2 = T2(-2.0) * ξ2
    d3 = ξ3 + T2(1.5)
    return N1, N2, N3, d1 * Δ1, d2 * Δ1, d3 * Δ1
end

@inline Base.@propagate_inbounds function topo_para(
    grid::    DeviceGrid{T1, T2}, 
    mpts::DeviceParticle{T1, T2},
    ix  ::Int
) where {T1, T2}
    # base index in the grid
    # note the base index needs a shift of 0.5h (see Taichi)
    mξx, mξy, mξz = mpts.ξ[ix, 1], mpts.ξ[ix, 2], mpts.ξ[ix, 3]
    bnx = unsafe_trunc(T1, floor((mξx - T2(0.5) * grid.h - grid.x1) * grid.invh)) # column
    bny = unsafe_trunc(T1, floor((mξy - T2(0.5) * grid.h - grid.y1) * grid.invh)) # row (自下而上)
    bnz = unsafe_trunc(T1, floor((mξz - T2(0.5) * grid.h - grid.z1) * grid.invh)) # layer
    bid = grid.nx * grid.ny * bnz + grid.ny * bnx + bny + 1
    gξx = grid.x1 + bnx * grid.h
    gξy = grid.y1 + bny * grid.h
    gξz = grid.z1 + bnz * grid.h
    # compute x value in the shape function N(x)
    x1 = grid.invh * (mξx - gξx)
    x2 = grid.invh * (mξx - gξx - grid.h)
    x3 = grid.invh * (mξx - gξx - grid.h - grid.h)
    y1 = grid.invh * (mξy - gξy)
    y2 = grid.invh * (mξy - gξy - grid.h)
    y3 = grid.invh * (mξy - gξy - grid.h - grid.h)
    z1 = grid.invh * (mξz - gξz)
    z2 = grid.invh * (mξz - gξz - grid.h)
    z3 = grid.invh * (mξz - gξz - grid.h - grid.h)
    # compute shape function value
    Nx1, Nx2, Nx3, dx1, dx2, dx3 = bspline2basis(x1, x2, x3, grid.invh)
    Ny1, Ny2, Ny3, dy1, dy2, dy3 = bspline2basis(y1, y2, y3, grid.invh)
    Nz1, Nz2, Nz3, dz1, dz2, dz3 = bspline2basis(z1, z2, z3, grid.invh)
    # assign the value (in order)
    Nxs, Nys, Nzs = (Nx1, Nx2, Nx3), (Ny1, Ny2, Ny3), (Nz1, Nz2, Nz3)
    dxs, dys, dzs = (dx1, dx2, dx3), (dy1, dy2, dy3), (dz1, dz2, dz3)
    return bid, Nxs, Nys, Nzs, dxs, dys, dzs
end

macro N∂N(grid, mpts, ix, body)
    g   = grid
    mp  = mpts
    ixv = ix
    b   = body

    return esc(quote
        bid, Nxs, Nys, Nzs, dxs, dys, dzs = MaterialPointSolver.topo_para($g, $mp, $ixv)

        @KAunroll for k = Int32(1):Int32(3)
                  for i = Int32(1):Int32(3)
                  for j = Int32(1):Int32(3)

                    p2n = Int32(bid +     (j - 1) +
                          $g.ny *         (i - 1) +
                          $g.ny * $g.nx * (k - 1))

                    Nij = Nxs[i] * Nys[j] * Nzs[k]
                    ∂Nx = dxs[i] * Nys[j] * Nzs[k]
                    ∂Ny = Nxs[i] * dys[j] * Nzs[k]
                    ∂Nz = Nxs[i] * Nys[j] * dzs[k]

                    $(b)
                end
            end
        end
    end)
end

macro Nij(grid, mpts, ix, body)
    g   = grid
    mp  = mpts
    ixv = ix
    b   = body

    return esc(quote
        bid, Nxs, Nys, Nzs, dxs, dys, dzs = MaterialPointSolver.topo_para($g, $mp, $ixv)

        @KAunroll for k = Int32(1):Int32(3)
                  for i = Int32(1):Int32(3)
                  for j = Int32(1):Int32(3)

                    p2n = Int32(bid +     (j - 1) +
                          $g.ny *         (i - 1) +
                          $g.ny * $g.nx * (k - 1))

                    Nij = Nxs[i] * Nys[j] * Nzs[k]

                    $(b)
                end
            end
        end
    end)
end