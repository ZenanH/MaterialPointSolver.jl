function resetgridstatus!(grid, ϵ::Precision{T1, T2}) where {T1, T2}
    fill!(grid.ms, T2(0.0))
    fill!(grid.vs, T2(0.0))
    fill!(grid.fs, T2(0.0))
    fill!(grid.ps, T2(0.0))
    fill!(grid.Ω , T2(0.0))
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

@kernel function P2G!(
    grid, mpts, gravity::T2, dim::Dim3D, ϵ::Precision{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        vol, ms = mpts.Ω[ix], mpts.Ω[ix] * mpts.ρs[ix]
        mppsx, mppsy, mppsz = mpts.vs[1, ix] * ms, mpts.vs[2, ix] * ms, mpts.vs[3, ix] * ms
        σxx, σyy, σzz = mpts.σij[1, ix], mpts.σij[2, ix], mpts.σij[3, ix]
        σxy, σyz, σzx = mpts.σij[4, ix], mpts.σij[5, ix], mpts.σij[6, ix]
        @KAunroll for iy in 1:mpts.NIC
            Nij = mpts.Nij[iy, ix]
            ∂Nx = mpts.∂Nx[iy, ix]
            ∂Ny = mpts.∂Ny[iy, ix]
            ∂Nz = mpts.∂Nz[iy, ix]
            p2n = mpts.p2n[iy, ix]
            # compute nodal mass
            @Σ grid.ms[p2n] += Nij * ms
            # compute nodal momentum
            @Σ grid.ps[1, p2n] += Nij * mppsx
            @Σ grid.ps[2, p2n] += Nij * mppsy
            @Σ grid.ps[3, p2n] += Nij * mppsz
            # compute nodal total force for solid
            @Σ grid.fs[1, p2n] += -vol * (∂Nx * σxx + ∂Ny * σxy + ∂Nz * σzx)
            @Σ grid.fs[2, p2n] += -vol * (∂Ny * σyy + ∂Nx * σxy + ∂Nz * σyz)
            @Σ grid.fs[3, p2n] += -vol * (∂Nz * σzz + ∂Nx * σzx + ∂Ny * σyz) + Nij * ms * gravity
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

@kernel function solvegrid!(
    grid, Δt::T2, dim::Dim3D, ϵ::Precision{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : 1 / grid.ms[ix]
        # boundary condition
        grid.vx_s_idx[ix] ≠ T1(0) ? grid.ps[1, ix] = grid.vx_s_val[ix] : nothing
        grid.vy_s_idx[ix] ≠ T1(0) ? grid.ps[2, ix] = grid.vy_s_val[ix] : nothing
        grid.vz_s_idx[ix] ≠ T1(0) ? grid.ps[3, ix] = grid.vz_s_val[ix] : nothing
        # compute nodal velocity
        grid.vs[1, ix] = grid.ps[1, ix] * ms_denom
        grid.vs[2, ix] = grid.ps[2, ix] * ms_denom
        grid.vs[3, ix] = grid.ps[3, ix] * ms_denom
        # update nodal velocity
        grid.vsT[1, ix] = grid.vs[1, ix] + grid.fs[1, ix] * Δt * ms_denom
        grid.vsT[2, ix] = grid.vs[2, ix] + grid.fs[2, ix] * Δt * ms_denom
        grid.vsT[3, ix] = grid.vs[3, ix] + grid.fs[3, ix] * Δt * ms_denom
        # boundary condition
        grid.vx_s_idx[ix] ≠ T1(0) ? grid.vsT[1, ix] = grid.vx_s_val[ix] : nothing
        grid.vy_s_idx[ix] ≠ T1(0) ? grid.vsT[2, ix] = grid.vy_s_val[ix] : nothing
        grid.vz_s_idx[ix] ≠ T1(0) ? grid.vsT[3, ix] = grid.vz_s_val[ix] : nothing
        # reset grid momentum
        grid.ps[1, ix] = T2(0.0)
        grid.ps[2, ix] = T2(0.0)
        grid.ps[3, ix] = T2(0.0)
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

@kernel function doublemapping1!(
    grid, mpts, Δt::T2, FLIP::T2, dim::Dim3D, ϵ::Precision{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        ξx = ξy = ξz = vx = vy = vz = T2(0.0)
        @KAunroll for iy in 1:mpts.NIC
            Nij = mpts.Nij[iy, ix]
            p2n = mpts.p2n[iy, ix]
            ξx += Nij *  grid.vsT[1, p2n]
            ξy += Nij *  grid.vsT[2, p2n]
            ξz += Nij *  grid.vsT[3, p2n]
            vx += Nij * (grid.vsT[1, p2n] - grid.vs[1, p2n])
            vy += Nij * (grid.vsT[2, p2n] - grid.vs[2, p2n])
            vz += Nij * (grid.vsT[3, p2n] - grid.vs[3, p2n])
        end
        # update particle position
        mpts.ξ[1, ix] += Δt * ξx
        mpts.ξ[2, ix] += Δt * ξy
        mpts.ξ[3, ix] += Δt * ξz
        # update particle velocity
        mpts.vs[1, ix] = FLIP * (mpts.vs[1, ix] + vx) + (T2(1.0) - FLIP) * ξx
        mpts.vs[2, ix] = FLIP * (mpts.vs[2, ix] + vy) + (T2(1.0) - FLIP) * ξy
        mpts.vs[3, ix] = FLIP * (mpts.vs[3, ix] + vz) + (T2(1.0) - FLIP) * ξz
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

@kernel function doublemapping2!(grid, mpts, dim::Dim3D, ϵ::Precision{T1, T2}) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        ms = mpts.ρs[ix] * mpts.Ω[ix]
        px, py, pz = ms * mpts.vs[1, ix], ms * mpts.vs[2, ix], ms * mpts.vs[3, ix]
        # update particle position & velocity
        for iy in Int32(1):Int32(mpts.NIC)
            Nij = mpts.Nij[iy, ix]
            p2n = mpts.p2n[iy, ix]
            @Σ grid.ps[1, p2n] += px * Nij
            @Σ grid.ps[2, p2n] += py * Nij
            @Σ grid.ps[3, p2n] += pz * Nij
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

@kernel function doublemapping3!(grid, dim::Dim3D, ϵ::Precision{T1, T2}) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        # compute nodal velocities
        grid.vs[1, ix] = grid.ps[1, ix] * ms_denom
        grid.vs[2, ix] = grid.ps[2, ix] * ms_denom
        grid.vs[3, ix] = grid.ps[3, ix] * ms_denom
        # fixed Dirichlet nodes
        grid.vx_s_idx[ix] ≠ T1(0) ? grid.vs[1, ix] = grid.vx_s_val[ix] : nothing
        grid.vy_s_idx[ix] ≠ T1(0) ? grid.vs[2, ix] = grid.vy_s_val[ix] : nothing
        grid.vz_s_idx[ix] ≠ T1(0) ? grid.vs[3, ix] = grid.vz_s_val[ix] : nothing
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

@kernel function G2P!(
    grid, mpts, Δt::T2, dim::Dim3D, ϵ::Precision{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        df1 = df2 = df3 = df4 = df5 = df6 = df7 = df8 = df9 = T2(0.0)
        @KAunroll for iy in 1:mpts.NIC
            p2n = mpts.p2n[iy, ix]
            ∂Nx = mpts.∂Nx[iy, ix]
            ∂Ny = mpts.∂Ny[iy, ix]
            ∂Nz = mpts.∂Nz[iy, ix]
            # compute solid incremental deformation gradient
            df1 += grid.vs[1, p2n] * ∂Nx
            df2 += grid.vs[1, p2n] * ∂Ny
            df3 += grid.vs[1, p2n] * ∂Nz
            df4 += grid.vs[2, p2n] * ∂Nx
            df5 += grid.vs[2, p2n] * ∂Ny
            df6 += grid.vs[2, p2n] * ∂Nz
            df7 += grid.vs[3, p2n] * ∂Nx
            df8 += grid.vs[3, p2n] * ∂Ny
            df9 += grid.vs[3, p2n] * ∂Nz
        end
        df1 *= Δt; df2 *= Δt; df3 *= Δt 
        df4 *= Δt; df5 *= Δt; df6 *= Δt
        df7 *= Δt; df8 *= Δt; df9 *= Δt
        # deformation gradient matrix
        update_F!(mpts, df1, df2, df3, df4, df5, df6, df7, df8, df9, T2(1.0), ix, dim)
        # update jacobian value and particle volume
        J = detF(mpts, ix, dim)
        mpts.Ω[ix] = J * mpts.Ω0[ix]
        mpts.ρs[ix] = mpts.ρs0[ix] / J
        liE!(mpts, df1, df2, df3, df4, df5, df6, df7, df8, df9, ix, dim, ϵ)
        dpP!(mpts, ix, dim, ϵ)
    end
end