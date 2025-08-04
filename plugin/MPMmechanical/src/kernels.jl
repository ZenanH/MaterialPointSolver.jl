function resetgridstatus!(grid, ϵ::Precision{T1, T2}) where {T1, T2}
    fill!(grid.ms, T2(0.0))
    fill!(grid.vs, T2(0.0))
    fill!(grid.fs, T2(0.0))
    fill!(grid.ps, T2(0.0))
    fill!(grid.Ω , T2(0.0))
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

@kernel function P2G!(
    grid, mpts, gravity::T2, dim::Dim2D, ϵ::Precision{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        vol, ms = mpts.Ω[ix], mpts.Ω[ix] * mpts.ρs[ix]
        mppsx, mppsy = mpts.vs[1, ix] * ms, mpts.vs[2, ix] * ms
        σxx, σyy, σxy = mpts.σij[1, ix], mpts.σij[2, ix], mpts.σij[4, ix]
        @KAunroll for iy in 1:mpts.NIC
            Nij = mpts.Nij[iy, ix]
            ∂Nx = mpts.∂Nx[iy, ix]
            ∂Ny = mpts.∂Ny[iy, ix]
            p2n = mpts.p2n[iy, ix]
            # compute nodal mass
            @Σ grid.ms[p2n] += Nij * ms
            # compute nodal momentum
            @Σ grid.ps[1, p2n] += Nij * mppsx
            @Σ grid.ps[2, p2n] += Nij * mppsy
            # compute nodal total force for solid
            @Σ grid.fs[1, p2n] += -vol * (∂Nx * σxx + ∂Ny * σxy)
            @Σ grid.fs[2, p2n] += -vol * (∂Ny * σyy + ∂Nx * σxy) + Nij * ms * gravity
        end
    end
end

@kernel function solvegrid!(
    grid, Δt::T2, dim::Dim2D, ϵ::Precision{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : 1 / grid.ms[ix]
        # boundary condition
        grid.vx_s_idx[ix] ≠ T1(0) ? grid.ps[1, ix] = grid.vx_s_val[ix] : nothing
        grid.vy_s_idx[ix] ≠ T1(0) ? grid.ps[2, ix] = grid.vy_s_val[ix] : nothing
        # compute nodal velocity
        grid.vs[1, ix] = grid.ps[1, ix] * ms_denom
        grid.vs[2, ix] = grid.ps[2, ix] * ms_denom
        # update nodal velocity
        grid.vsT[1, ix] = grid.vs[1, ix] + grid.fs[1, ix] * Δt * ms_denom
        grid.vsT[2, ix] = grid.vs[2, ix] + grid.fs[2, ix] * Δt * ms_denom
        # boundary condition
        grid.vx_s_idx[ix] ≠ T1(0) ? grid.vsT[1, ix] = grid.vx_s_val[ix] : nothing
        grid.vy_s_idx[ix] ≠ T1(0) ? grid.vsT[2, ix] = grid.vy_s_val[ix] : nothing
        # reset grid momentum
        grid.ps[1, ix] = T2(0.0)
        grid.ps[2, ix] = T2(0.0)
    end
end

@kernel function doublemapping1!(
    grid, mpts, Δt::T2, FLIP::T2, dim::Dim2D, ϵ::Precision{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        ξx = ξy = vx = vy = T2(0.0)
        @KAunroll for iy in 1:mpts.NIC
            Nij = mpts.Nij[iy, ix]
            p2n = mpts.p2n[iy, ix]
            ξx += Nij *  grid.vsT[1, p2n]
            ξy += Nij *  grid.vsT[2, p2n]
            vx += Nij * (grid.vsT[1, p2n] - grid.vs[1, p2n])
            vy += Nij * (grid.vsT[2, p2n] - grid.vs[2, p2n])
        end
        # update particle position
        mpts.ξ[1, ix] += Δt * ξx
        mpts.ξ[2, ix] += Δt * ξy
        # update particle velocity
        mpts.vs[1, ix] = FLIP * (mpts.vs[1, ix] + vx) + (T2(1.0) - FLIP) * ξx
        mpts.vs[2, ix] = FLIP * (mpts.vs[2, ix] + vy) + (T2(1.0) - FLIP) * ξy
    end
end

@kernel function doublemapping2!(grid, mpts, dim::Dim2D, ϵ::Precision{T1, T2}) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        ms = mpts.ρs[ix] * mpts.Ω[ix]
        px, py = ms * mpts.vs[1, ix], ms * mpts.vs[2, ix]
        # update particle position & velocity
        for iy in Int32(1):Int32(mpts.NIC)
            Nij = mpts.Nij[iy, ix]
            p2n = mpts.p2n[iy, ix]
            @Σ grid.ps[1, p2n] += px * Nij
            @Σ grid.ps[2, p2n] += py * Nij
        end
    end
end

@kernel function doublemapping3!(grid, dim::Dim2D, ϵ::Precision{T1, T2}) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        # compute nodal velocities
        grid.vs[1, ix] = grid.ps[1, ix] * ms_denom
        grid.vs[2, ix] = grid.ps[2, ix] * ms_denom
        # fixed Dirichlet nodes
        grid.vx_s_idx[ix] ≠ T1(0) ? grid.vs[1, ix] = grid.vx_s_val[ix] : nothing
        grid.vy_s_idx[ix] ≠ T1(0) ? grid.vs[2, ix] = grid.vy_s_val[ix] : nothing
    end
end

@kernel function G2P!(
    grid, mpts, Δt::T2, dim::Dim2D, ϵ::Precision{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        df1 = df2 = df3 = df4 = T2(0.0)
        @KAunroll for iy in 1:mpts.NIC
            p2n = mpts.p2n[iy, ix]
            ∂Nx = mpts.∂Nx[iy, ix]
            ∂Ny = mpts.∂Ny[iy, ix]
            # compute solid incremental deformation gradient
            df1 += grid.vs[1, p2n] * ∂Nx
            df2 += grid.vs[1, p2n] * ∂Ny
            df3 += grid.vs[2, p2n] * ∂Nx
            df4 += grid.vs[2, p2n] * ∂Ny
        end
        df1 *= Δt; df2 *= Δt; df3 *= Δt; df4 *= Δt
        # deformation gradient matrix
        update_F!(mpts, df1, df2, df3, df4, T2(1.0), ix, dim)
        # update jacobian value and particle volume
        J = detF(mpts, ix, dim)
        mpts.Ω[ix] = J * mpts.Ω0[ix]
        mpts.ρs[ix] = mpts.ρs0[ix] / J
        liE!(mpts, df1, df2, df3, df4, ix, dim, ϵ)
        dpP!(mpts, ix, dim, ϵ)
    end
end

function launchkernels!(conf, dev_grid, dev_mpts, Δt)
    dev   = conf.dev
    basis = conf.basis
    dim   = conf.dim
    ϵ     = conf.ϵ
    G     = dev_mpts.G
    FLIP  = dev_mpts.FLIP
    resetgridstatus!(dev_grid, ϵ)
    resetmpstatus!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, basis, dim, ϵ)
    P2G!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, G, dim, ϵ)
    solvegrid!(dev)(ndrange=dev_grid.ni, dev_grid, Δt, dim, ϵ)
    doublemapping1!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, Δt, FLIP, dim, ϵ)
    doublemapping2!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, dim, ϵ)
    doublemapping3!(dev)(ndrange=dev_grid.ni, dev_grid, dim, ϵ)
    G2P!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, Δt, dim, ϵ)
end