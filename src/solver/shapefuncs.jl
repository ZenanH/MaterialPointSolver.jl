export @N∂NBspline2, @NijBspline2
export @N∂NuGIMP, @NijuGIMP

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

@inline Base.@propagate_inbounds function uGIMPbasis(
    x1::T2, 
    x2::T2, 
    x3::T2, 
    h ::T2,
    lp::T2
) where T2
    N1 = d1 = N2 = d2 = N3 = d3 = T2(0.0)
    absx1, absx2, absx3 = abs(x1), abs(x2), abs(x3)

    if absx1 < T2(0.5)*lp
        N1 = T2(1.0) - ((T2(4.0) * x1 * x1 + lp * lp) / (T2(4.0) * h * lp))
        d1 = -((T2(8.0) * x1) / (T2(4.0) * h * lp))
    elseif (T2(0.5) * lp) ≤ absx1 < (h - T2(0.5) * lp)
        N1 = T2(1.0) - (absx1 / h)
        d1 = sign(x1) * (T2(-1.0) / h)
    elseif (h - T2(0.5) * lp) ≤ absx1 < (h + T2(0.5) * lp)
        N1 = ((h + T2(0.5) * lp - absx1) * (h + T2(0.5) * lp - absx1)) / (T2(2.0) * h * lp)
        d1 = -sign(x1) * ((h + T2(0.5) * lp - absx1) / (h * lp))
    end

    if absx2 < T2(0.5)*lp
        N2 = T2(1.0) - ((T2(4.0) * x2 * x2 + lp * lp) / (T2(4.0) * h * lp))
        d2 = -((T2(8.0) * x2) / (T2(4.0) * h * lp))
    elseif (T2(0.5) * lp) ≤ absx2 < (h - T2(0.5) * lp)
        N2 = T2(1.0) - (absx2 / h)
        d2 = sign(x2) * (T2(-1.0) / h)
    elseif (h - T2(0.5) * lp) ≤ absx2 < (h + T2(0.5) * lp)
        N2 = ((h + T2(0.5) * lp - absx2) * (h + T2(0.5) * lp - absx2)) / (T2(2.0) * h * lp)
        d2 = -sign(x2) * ((h + T2(0.5) * lp - absx2) / (h * lp))
    end

    if absx3 < T2(0.5)*lp
        N3 = T2(1.0) - ((T2(4.0) * x3 * x3 + lp * lp) / (T2(4.0) * h * lp))
        d3 = -((T2(8.0) * x3) / (T2(4.0) * h * lp))
    elseif (T2(0.5) * lp) ≤ absx3 < (h - T2(0.5) * lp)
        N3 = T2(1.0) - (absx3 / h)
        d3 = sign(x3) * (T2(-1.0) / h)
    elseif (h - T2(0.5) * lp) ≤ absx3 < (h + T2(0.5) * lp)
        N3 = ((h + T2(0.5) * lp - absx3) * (h + T2(0.5) * lp - absx3)) / (T2(2.0) * h * lp)
        d3 = -sign(x3) * ((h + T2(0.5) * lp - absx3) / (h * lp))
    end
    
    return N1, N2, N3, d1, d2, d3
end

@inline Base.@propagate_inbounds function topo_para_Bspline2(
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

@inline Base.@propagate_inbounds function topo_para_uGIMP(
    grid::    DeviceGrid{T1, T2}, 
    mpts::DeviceParticle{T1, T2},
    ix  ::Int
) where {T1, T2}
    # base index in the grid
    # note the base index needs a shift of 0.5h (see Taichi)
    mξx, mξy, mξz = mpts.ξ[ix, 1], mpts.ξ[ix, 2], mpts.ξ[ix, 3]
    bnx = unsafe_trunc(T1, floor((mξx - T2(0.5) * mpts.h - grid.x1) * grid.invh))
    bny = unsafe_trunc(T1, floor((mξy - T2(0.5) * mpts.h - grid.y1) * grid.invh))
    bnz = unsafe_trunc(T1, floor((mξz - T2(0.5) * mpts.h - grid.z1) * grid.invh))
    bid = grid.nx * grid.ny * bnz + grid.ny * bnx + bny + 1
    gξx = grid.x1 + bnx * grid.h
    gξy = grid.y1 + bny * grid.h
    gξz = grid.z1 + bnz * grid.h
    # compute x value in the basis function N(x)
    x1 = mξx - gξx
    x2 = mξx - gξx - grid.h
    x3 = mξx - gξx - grid.h - grid.h
    y1 = mξy - gξy
    y2 = mξy - gξy - grid.h
    y3 = mξy - gξy - grid.h - grid.h
    z1 = mξz - gξz
    z2 = mξz - gξz - grid.h
    z3 = mξz - gξz - grid.h - grid.h
    Nx1, Nx2, Nx3, dx1, dx2, dx3 = uGIMPbasis(x1, x2, x3, grid.h, mpts.h)
    Ny1, Ny2, Ny3, dy1, dy2, dy3 = uGIMPbasis(y1, y2, y3, grid.h, mpts.h)
    Nz1, Nz2, Nz3, dz1, dz2, dz3 = uGIMPbasis(z1, z2, z3, grid.h, mpts.h)
    # assign the value (in order)
    Nxs, Nys, Nzs = (Nx1, Nx2, Nx3), (Ny1, Ny2, Ny3), (Nz1, Nz2, Nz3)
    dxs, dys, dzs = (dx1, dx2, dx3), (dy1, dy2, dy3), (dz1, dz2, dz3)
    return bid, Nxs, Nys, Nzs, dxs, dys, dzs
end

macro N∂NBspline2(grid, mpts, ix, body)
    g   = grid
    mp  = mpts
    ixv = ix
    b   = body

    return esc(quote
        bid, Nxs, Nys, Nzs, dxs, dys, dzs = MaterialPointSolver.topo_para_Bspline2($g, $mp, $ixv)

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

macro N∂NuGIMP(grid, mpts, ix, body)
    g   = grid
    mp  = mpts
    ixv = ix
    b   = body

    return esc(quote
        bid, Nxs, Nys, Nzs, dxs, dys, dzs = MaterialPointSolver.topo_para_uGIMP($g, $mp, $ixv)

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

macro NijBspline2(grid, mpts, ix, body)
    g   = grid
    mp  = mpts
    ixv = ix
    b   = body

    return esc(quote
        bid, Nxs, Nys, Nzs, dxs, dys, dzs = MaterialPointSolver.topo_para_Bspline2($g, $mp, $ixv)

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

macro NijuGIMP(grid, mpts, ix, body)
    g   = grid
    mp  = mpts
    ixv = ix
    b   = body

    return esc(quote
        bid, Nxs, Nys, Nzs, dxs, dys, dzs = MaterialPointSolver.topo_para_uGIMP($g, $mp, $ixv)

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