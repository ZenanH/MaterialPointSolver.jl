export resetgridstatus!
export resetmpstatus!
export p2g!
export solvegrid!
export doublemapping1!
export doublemapping2!
export doublemapping3!
export g2p!

function resetgridstatus!(grid::DeviceGrid{T1, T2}) where {T1, T2}
    fill!(grid.ms, T2(0.0))
    fill!(grid.fs, T2(0.0))
    fill!(grid.ps, T2(0.0))
    fill!(grid.vs, T2(0.0))
end

@kernel function resetmpstatus!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}, ::Bspline2Basis) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
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
        # compute x value in the basis function N(x)
        x1 = grid.invh * (mξx - gξx)
        x2 = grid.invh * (mξx - gξx - grid.h)
        x3 = grid.invh * (mξx - gξx - grid.h - grid.h)
        y1 = grid.invh * (mξy - gξy)
        y2 = grid.invh * (mξy - gξy - grid.h)
        y3 = grid.invh * (mξy - gξy - grid.h - grid.h)
        z1 = grid.invh * (mξz - gξz)
        z2 = grid.invh * (mξz - gξz - grid.h)
        z3 = grid.invh * (mξz - gξz - grid.h - grid.h)
        Nx1, Nx2, Nx3, dx1, dx2, dx3 = bspline2basis(x1, x2, x3, grid.invh)
        Ny1, Ny2, Ny3, dy1, dy2, dy3 = bspline2basis(y1, y2, y3, grid.invh)
        Nz1, Nz2, Nz3, dz1, dz2, dz3 = bspline2basis(z1, z2, z3, grid.invh)
        # assign the value (in order)
        it = Int32(1)
        Nxs, Nys, Nzs = (Nx1, Nx2, Nx3), (Ny1, Ny2, Ny3), (Nz1, Nz2, Nz3)
        dxs, dys, dzs = (dx1, dx2, dx3), (dy1, dy2, dy3), (dz1, dz2, dz3)
        @KAunroll for k in Int32(1):Int32(3) # z-direction
            for i in Int32(1):Int32(3)       # x-direction
                for j in Int32(1):Int32(3)   # y-direction
                    mpts.p2n[ix, it] = bid +               (j - Int32(1)) + 
                                       grid.ny *           (i - Int32(1)) + 
                                       grid.ny * grid.nx * (k - Int32(1))
                    mpts.Nij[ix, it] = Nxs[i] * Nys[j] * Nzs[k]
                    mpts.∂Nx[ix, it] = dxs[i] * Nys[j] * Nzs[k]
                    mpts.∂Ny[ix, it] = Nxs[i] * dys[j] * Nzs[k]
                    mpts.∂Nz[ix, it] = Nxs[i] * Nys[j] * dzs[k]
                    it += Int32(1)
                end
            end
        end
    end
end

@kernel function resetmpstatus!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}, ::LinearBasis) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        # base index in the grid
        # note the base index needs a shift of 0.5h (see Taichi)
        mξx, mξy, mξz = mpts.ξ[ix, 1], mpts.ξ[ix, 2], mpts.ξ[ix, 3]
        bnx = unsafe_trunc(T1, floor((mξx - grid.x1) * grid.invh))
        bny = unsafe_trunc(T1, floor((mξy - grid.y1) * grid.invh))
        bnz = unsafe_trunc(T1, floor((mξz - grid.z1) * grid.invh))
        bid = grid.nx * grid.ny * bnz + grid.ny * bnx + bny + 1
        gξx = grid.x1 + bnx * grid.h
        gξy = grid.y1 + bny * grid.h
        gξz = grid.z1 + bnz * grid.h
        # compute x value in the basis function N(x)
        x1, y1, z1 = mξx - gξx, mξy - gξy, mξz - gξz
        x2, y2, z2 = grid.h - x1, grid.h - y1, grid.h - z1
        Nx1, Nx2, ∂x1, ∂x2 = linearbasis(x1, x2, grid.invh)
        Ny1, Ny2, ∂y1, ∂y2 = linearbasis(y1, y2, grid.invh)
        Nz1, Nz2, ∂z1, ∂z2 = linearbasis(z1, z2, grid.invh)
        # assign the value (in order)
        it = Int32(1)
        Nxs, Nys, Nzs = (Nx1, Nx2), (Ny1, Ny2), (Nz1, Nz2)
        dxs, dys, dzs = (∂x1, ∂x2), (∂y1, ∂y2), (∂z1, ∂z2)
        @KAunroll for k in Int32(1):Int32(2) # z-direction
            for i in Int32(1):Int32(2)       # x-direction
                for j in Int32(1):Int32(2)   # y-direction
                    mpts.p2n[ix, it] = bid +               (j - Int32(1)) + 
                                       grid.ny *           (i - Int32(1)) + 
                                       grid.ny * grid.nx * (k - Int32(1))
                    mpts.Nij[ix, it] = Nxs[i] * Nys[j] * Nzs[k]
                    mpts.∂Nx[ix, it] = dxs[i] * Nys[j] * Nzs[k]
                    mpts.∂Ny[ix, it] = Nxs[i] * dys[j] * Nzs[k]
                    mpts.∂Nz[ix, it] = Nxs[i] * Nys[j] * dzs[k]
                    it += Int32(1)
                end
            end
        end
    end
end

@kernel function resetmpstatus!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}, ::uGIMPBasis) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
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
        it = Int32(1)
        Nxs, Nys, Nzs = (Nx1, Nx2, Nx3), (Ny1, Ny2, Ny3), (Nz1, Nz2, Nz3)
        dxs, dys, dzs = (dx1, dx2, dx3), (dy1, dy2, dy3), (dz1, dz2, dz3)
        @KAunroll for k in Int32(1):Int32(3) # z-direction
            for i in Int32(1):Int32(3)       # x-direction
                for j in Int32(1):Int32(3)   # y-direction
                    mpts.p2n[ix, it] = bid +               (j - Int32(1)) + 
                                       grid.ny *           (i - Int32(1)) + 
                                       grid.ny * grid.nx * (k - Int32(1))
                    mpts.Nij[ix, it] = Nxs[i] * Nys[j] * Nzs[k]
                    mpts.∂Nx[ix, it] = dxs[i] * Nys[j] * Nzs[k]
                    mpts.∂Ny[ix, it] = Nxs[i] * dys[j] * Nzs[k]
                    mpts.∂Nz[ix, it] = Nxs[i] * Nys[j] * dzs[k]
                    it += Int32(1)
                end
            end
        end
    end
end

@kernel function resetmpstatus!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}, ::Bspline3Basis) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        # base index in the grid
        # note the base index needs a shift of 1h (see Taichi)
        mξx, mξy, mξz = mpts.ξ[ix, 1], mpts.ξ[ix, 2], mpts.ξ[ix, 3]
        bnx = unsafe_trunc(T1, floor((mξx - grid.h - grid.x1) * grid.invh))
        bny = unsafe_trunc(T1, floor((mξy - grid.h - grid.y1) * grid.invh))
        bnz = unsafe_trunc(T1, floor((mξz - grid.h - grid.z1) * grid.invh))
        bid = grid.nx * grid.ny * bnz + grid.ny * bnx + bny + 1
        gξx = grid.x1 + bnx * grid.h
        gξy = grid.y1 + bny * grid.h
        gξz = grid.z1 + bnz * grid.h
        # compute x value in the basis function N(x)
        rx1 = (mξx - gξx                           ) * grid.invh
        rx2 = (mξx - gξx - grid.h                  ) * grid.invh
        rx3 = (mξx - gξx - grid.h - grid.h         ) * grid.invh
        rx4 = (mξx - gξx - grid.h - grid.h - grid.h) * grid.invh
        ry1 = (mξy - gξy                           ) * grid.invh
        ry2 = (mξy - gξy - grid.h                  ) * grid.invh
        ry3 = (mξy - gξy - grid.h - grid.h         ) * grid.invh
        ry4 = (mξy - gξy - grid.h - grid.h - grid.h) * grid.invh
        rz1 = (mξz - gξz                           ) * grid.invh
        rz2 = (mξz - gξz - grid.h                  ) * grid.invh
        rz3 = (mξz - gξz - grid.h - grid.h         ) * grid.invh
        rz4 = (mξz - gξz - grid.h - grid.h - grid.h) * grid.invh
        tx1 = get_type(gξx                           , grid.x1, grid.x2, grid.h)
        tx2 = get_type(gξx + grid.h                  , grid.x1, grid.x2, grid.h)
        tx3 = get_type(gξx + grid.h + grid.h         , grid.x1, grid.x2, grid.h)
        tx4 = get_type(gξx + grid.h + grid.h + grid.h, grid.x1, grid.x2, grid.h)
        ty1 = get_type(gξy                           , grid.y1, grid.y2, grid.h)
        ty2 = get_type(gξy + grid.h                  , grid.y1, grid.y2, grid.h)
        ty3 = get_type(gξy + grid.h + grid.h         , grid.y1, grid.y2, grid.h)
        ty4 = get_type(gξy + grid.h + grid.h + grid.h, grid.y1, grid.y2, grid.h)
        tz1 = get_type(gξz                           , grid.z1, grid.z2, grid.h)
        tz2 = get_type(gξz + grid.h                  , grid.z1, grid.z2, grid.h)
        tz3 = get_type(gξz + grid.h + grid.h         , grid.z1, grid.z2, grid.h)
        tz4 = get_type(gξz + grid.h + grid.h + grid.h, grid.z1, grid.z2, grid.h)
        Nx1, ∂x1 = bspline3basis(rx1, grid.h, tx1)
        Nx2, ∂x2 = bspline3basis(rx2, grid.h, tx2)
        Nx3, ∂x3 = bspline3basis(rx3, grid.h, tx3)
        Nx4, ∂x4 = bspline3basis(rx4, grid.h, tx4)
        Ny1, ∂y1 = bspline3basis(ry1, grid.h, ty1)
        Ny2, ∂y2 = bspline3basis(ry2, grid.h, ty2)
        Ny3, ∂y3 = bspline3basis(ry3, grid.h, ty3)
        Ny4, ∂y4 = bspline3basis(ry4, grid.h, ty4)
        Nz1, ∂z1 = bspline3basis(rz1, grid.h, tz1)
        Nz2, ∂z2 = bspline3basis(rz2, grid.h, tz2)
        Nz3, ∂z3 = bspline3basis(rz3, grid.h, tz3)
        Nz4, ∂z4 = bspline3basis(rz4, grid.h, tz4)
        # assign the value (in order)
        it = Int32(1)
        Nxs, Nys, Nzs = (Nx1, Nx2, Nx3, Nx4), (Ny1, Ny2, Ny3, Ny4), (Nz1, Nz2, Nz3, Nz4)
        dxs, dys, dzs = (∂x1, ∂x2, ∂x3, ∂x4), (∂y1, ∂y2, ∂y3, ∂y4), (∂z1, ∂z2, ∂z3, ∂z4)
        @KAunroll for k in Int32(1):Int32(4) # z-direction
            for i in Int32(1):Int32(4)       # x-direction
                for j in Int32(1):Int32(4)   # y-direction
                    mpts.p2n[ix, it] = bid +               (j - Int32(1)) + 
                                       grid.ny *           (i - Int32(1)) + 
                                       grid.ny * grid.nx * (k - Int32(1))
                    mpts.Nij[ix, it] = Nxs[i] * Nys[j] * Nzs[k]
                    mpts.∂Nx[ix, it] = dxs[i] * Nys[j] * Nzs[k]
                    mpts.∂Ny[ix, it] = Nxs[i] * dys[j] * Nzs[k]
                    mpts.∂Nz[ix, it] = Nxs[i] * Nys[j] * dzs[k]
                    it += Int32(1)
                end
            end
        end
    end
end

@kernel function p2g!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        Ω, ms = mpts.Ω[ix], mpts.ρs[ix] * mpts.Ω[ix]
        msG = ms * mpts.G
        psx, psy, psz = mpts.vs[ix, 1] * ms, mpts.vs[ix, 2] * ms, mpts.vs[ix, 3] * ms
        σxx, σyy, σzz = mpts.σij[ix, 1], mpts.σij[ix, 2], mpts.σij[ix, 3]
        σxy, σyz, σzx = mpts.σij[ix, 4], mpts.σij[ix, 5], mpts.σij[ix, 6]
        @KAunroll for iy in 1:mpts.NIC
            Nij = mpts.Nij[ix, iy]
            ∂Nx = mpts.∂Nx[ix, iy]
            ∂Ny = mpts.∂Ny[ix, iy]
            ∂Nz = mpts.∂Nz[ix, iy]
            p2n = mpts.p2n[ix, iy]
            # compute nodal mass
            @Σ grid.ms[p2n] += Nij * ms
            # compute nodal momentum
            @Σ grid.vs[p2n, 1] += Nij * psx
            @Σ grid.vs[p2n, 2] += Nij * psy
            @Σ grid.vs[p2n, 3] += Nij * psz
            # compute nodal total force for solid
            @Σ grid.fs[p2n, 1] += -Ω * (∂Nx * σxx + ∂Ny * σxy + ∂Nz * σzx)
            @Σ grid.fs[p2n, 2] += -Ω * (∂Nx * σxy + ∂Ny * σyy + ∂Nz * σyz)
            @Σ grid.fs[p2n, 3] += -Ω * (∂Nx * σzx + ∂Ny * σyz + ∂Nz * σzz) + Nij * msG
        end
    end
end

@kernel function solvegrid!(grid::DeviceGrid{T1, T2}, Δt::T2) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : 1 / grid.ms[ix]
        # boundary condition
        grid.vsxi[ix] ≠ T1(0) ? grid.vs[ix, 1] = grid.vsxv[ix] : nothing
        grid.vsyi[ix] ≠ T1(0) ? grid.vs[ix, 2] = grid.vsyv[ix] : nothing
        grid.vszi[ix] ≠ T1(0) ? grid.vs[ix, 3] = grid.vszv[ix] : nothing
        # compute nodal velocity
        grid.vs[ix, 1] *= ms_denom
        grid.vs[ix, 2] *= ms_denom
        grid.vs[ix, 3] *= ms_denom
        # update nodal velocity
        grid.vsT[ix, 1] = grid.vs[ix, 1] + grid.fs[ix, 1] * Δt * ms_denom
        grid.vsT[ix, 2] = grid.vs[ix, 2] + grid.fs[ix, 2] * Δt * ms_denom
        grid.vsT[ix, 3] = grid.vs[ix, 3] + grid.fs[ix, 3] * Δt * ms_denom
        # boundary condition
        grid.vsxi[ix] ≠ T1(0) ? grid.vsT[ix, 1] = grid.vsxv[ix] : nothing
        grid.vsyi[ix] ≠ T1(0) ? grid.vsT[ix, 2] = grid.vsyv[ix] : nothing
        grid.vszi[ix] ≠ T1(0) ? grid.vsT[ix, 3] = grid.vszv[ix] : nothing
    end
end

@kernel function doublemapping1!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}, Δt::T2) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        FLIP = mpts.FLIP
        ξx = ξy = ξz = vx = vy = vz = T2(0.0)
        @KAunroll for iy in 1:mpts.NIC
            Nij = mpts.Nij[ix, iy]
            p2n = mpts.p2n[ix, iy]
            ξx += Nij *  grid.vsT[p2n, 1]
            ξy += Nij *  grid.vsT[p2n, 2]
            ξz += Nij *  grid.vsT[p2n, 3]
            vx += Nij * (grid.vsT[p2n, 1] - grid.vs[p2n, 1])
            vy += Nij * (grid.vsT[p2n, 2] - grid.vs[p2n, 2])
            vz += Nij * (grid.vsT[p2n, 3] - grid.vs[p2n, 3])
        end
        # update particle position
        mpts.ξ[ix, 1] += Δt * ξx
        mpts.ξ[ix, 2] += Δt * ξy
        mpts.ξ[ix, 3] += Δt * ξz
        # update particle velocity
        mpts.vs[ix, 1] = FLIP * (mpts.vs[ix, 1] + vx) + (T2(1.0) - FLIP) * ξx
        mpts.vs[ix, 2] = FLIP * (mpts.vs[ix, 2] + vy) + (T2(1.0) - FLIP) * ξy
        mpts.vs[ix, 3] = FLIP * (mpts.vs[ix, 3] + vz) + (T2(1.0) - FLIP) * ξz
    end
end

@kernel function doublemapping2!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        ms = mpts.ρs[ix] * mpts.Ω[ix]
        px, py, pz = ms * mpts.vs[ix, 1], ms * mpts.vs[ix, 2], ms * mpts.vs[ix, 3]
        # update particle position & velocity
        @KAunroll for iy in 1:mpts.NIC
            Nij = mpts.Nij[ix, iy]
            p2n = mpts.p2n[ix, iy]
            @Σ grid.ps[p2n, 1] += px * Nij
            @Σ grid.ps[p2n, 2] += py * Nij
            @Σ grid.ps[p2n, 3] += pz * Nij
        end
    end
end

@kernel function doublemapping3!(grid::DeviceGrid{T1, T2}, Δt::T2) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : 1 / grid.ms[ix]
        # compute nodal velocities
        grid.ps[ix, 1] = grid.ps[ix, 1] * ms_denom * Δt
        grid.ps[ix, 2] = grid.ps[ix, 2] * ms_denom * Δt
        grid.ps[ix, 3] = grid.ps[ix, 3] * ms_denom * Δt
        # fixed Dirichlet nodes
        grid.vsxi[ix] ≠ T1(0) ? grid.ps[ix, 1] = grid.vsxv[ix] : nothing
        grid.vsyi[ix] ≠ T1(0) ? grid.ps[ix, 2] = grid.vsyv[ix] : nothing
        grid.vszi[ix] ≠ T1(0) ? grid.ps[ix, 3] = grid.vszv[ix] : nothing
    end
end

@kernel function g2p!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        df1 = df2 = df3 = df4 = df5 = df6 = df7 = df8 = df9 = T2(0.0)
        @KAunroll for iy in 1:mpts.NIC
            p2n = mpts.p2n[ix, iy]
            ∂Nx = mpts.∂Nx[ix, iy]; ux = grid.ps[p2n, 1]
            ∂Ny = mpts.∂Ny[ix, iy]; uy = grid.ps[p2n, 2]
            ∂Nz = mpts.∂Nz[ix, iy]; uz = grid.ps[p2n, 3]
            # compute solid incremental deformation gradient
            df1 += ux * ∂Nx; df4 += uy * ∂Nx; df7 += uz * ∂Nx
            df2 += ux * ∂Ny; df5 += uy * ∂Ny; df8 += uz * ∂Ny
            df3 += ux * ∂Nz; df6 += uy * ∂Nz; df9 += uz * ∂Nz
        end
        # deformation gradient matrix
        update_F!(mpts, df1, df2, df3, df4, df5, df6, df7, df8, df9, ix)
        # update jacobian value and particle volume
        J = detF(mpts, ix)
        mpts.Ω[ix] = J * mpts.Ω0[ix]
        mpts.ρs[ix] = mpts.ρs0[ix] / J
        liE!(mpts, df1, df2, df3, df4, df5, df6, df7, df8, df9, ix)
        dpP!(mpts, ix)
    end
end