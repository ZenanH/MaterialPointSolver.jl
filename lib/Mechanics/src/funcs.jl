function resetgridstatus!(grid::DeviceGrid{T1, T2}) where {T1, T2}
    fill!(grid.ms, T2(0.0))
    fill!(grid.fs, T2(0.0))
    fill!(grid.ps, T2(0.0))
    fill!(grid.vs, T2(0.0))
end

@kernel function p2g!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        Ω, ms = mpts.Ω[ix], mpts.ρs[ix] * mpts.Ω[ix]
        msG = ms * mpts.G
        psx, psy, psz = mpts.vs[ix, 1] * ms, mpts.vs[ix, 2] * ms, mpts.vs[ix, 3] * ms
        σxx, σyy, σzz = mpts.σij[ix, 1], mpts.σij[ix, 2], mpts.σij[ix, 3]
        σxy, σyz, σzx = mpts.σij[ix, 4], mpts.σij[ix, 5], mpts.σij[ix, 6]
        @N∂NBspline2(grid, mpts, ix, begin
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
        end)
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
        # damping force for solid
        dampvs = -grid.ζs * sqrt(grid.fs[ix, 1] * grid.fs[ix, 1] + 
                                 grid.fs[ix, 2] * grid.fs[ix, 2] + 
                                 grid.fs[ix, 3] * grid.fs[ix, 3])
        # compute nodal total force for mixture
        Fs_x = grid.fs[ix, 1] + dampvs * sign(grid.vs[ix, 1])
        Fs_y = grid.fs[ix, 2] + dampvs * sign(grid.vs[ix, 2])
        Fs_z = grid.fs[ix, 3] + dampvs * sign(grid.vs[ix, 3])
        # update nodal velocity
        grid.vsT[ix, 1] = grid.vs[ix, 1] + Fs_x * Δt * ms_denom
        grid.vsT[ix, 2] = grid.vs[ix, 2] + Fs_y * Δt * ms_denom
        grid.vsT[ix, 3] = grid.vs[ix, 3] + Fs_z * Δt * ms_denom
        # boundary condition
        grid.vsxi[ix] ≠ T1(0) ? grid.vsT[ix, 1] = grid.vsxv[ix] : nothing
        grid.vsyi[ix] ≠ T1(0) ? grid.vsT[ix, 2] = grid.vsyv[ix] : nothing
        grid.vszi[ix] ≠ T1(0) ? grid.vsT[ix, 3] = grid.vszv[ix] : nothing
        # compute velocity increment
        grid.vs[ix, 1] = grid.vsT[ix, 1] - grid.vs[ix, 1]
        grid.vs[ix, 2] = grid.vsT[ix, 2] - grid.vs[ix, 2]
        grid.vs[ix, 3] = grid.vsT[ix, 3] - grid.vs[ix, 3]
    end
end

@kernel function doublemapping1!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        FLIP = mpts.FLIP
        ξx = ξy = ξz = vx = vy = vz = T2(0.0)
        @NijBspline2(grid, mpts, ix, begin
            ξx += Nij * grid.vsT[p2n, 1]
            ξy += Nij * grid.vsT[p2n, 2]
            ξz += Nij * grid.vsT[p2n, 3]
            vx += Nij * grid.vs[p2n, 1]
            vy += Nij * grid.vs[p2n, 2]
            vz += Nij * grid.vs[p2n, 3]
        end)
        # update particle velocity
        mpts.vs[ix, 1] = FLIP * (mpts.vs[ix, 1] + vx) + (T2(1.0) - FLIP) * ξx
        mpts.vs[ix, 2] = FLIP * (mpts.vs[ix, 2] + vy) + (T2(1.0) - FLIP) * ξy
        mpts.vs[ix, 3] = FLIP * (mpts.vs[ix, 3] + vz) + (T2(1.0) - FLIP) * ξz
        # update grid momentum
        ms = mpts.ρs[ix] * mpts.Ω[ix]
        px, py, pz = ms * mpts.vs[ix, 1], ms * mpts.vs[ix, 2], ms * mpts.vs[ix, 3]
        # update CFL conditions
        nid  = mpts.nid[ix]
        Ks   = mpts.Ks[nid]
        Gs   = mpts.Gs[nid]
        cdil = sqrt((Ks + T2(1.333333) * Gs) / mpts.ρs[ix]) # 4/3 ≈ 1.333333
        mpts.cfl[ix] = min(grid.h / (cdil + abs(mpts.vs[ix, 1])), 
                           grid.h / (cdil + abs(mpts.vs[ix, 2]))) 
        # update grid momentum
        @NijBspline2(grid, mpts, ix, begin
            @Σ grid.ps[p2n, 1] += px * Nij
            @Σ grid.ps[p2n, 2] += py * Nij
            @Σ grid.ps[p2n, 3] += pz * Nij
        end)
    end
end

@kernel function doublemapping2!(grid::DeviceGrid{T1, T2}, Δt::T2) where {T1, T2}
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

@kernel function g2p!(
    grid ::    DeviceGrid{T1, T2}, 
    mpts ::DeviceParticle{T1, T2}, 
    t_eld::T2, 
    t_cur::T2, 
    Δt   ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        df1 = df2 = df3 = df4 = df5 = df6 = df7 = df8 = df9 = T2(0.0)
        ξx = ξy = ξz = T2(0.0)
        @N∂NBspline2(grid, mpts, ix, begin
            ux = grid.ps[p2n, 1]; ξx += ux * Nij
            uy = grid.ps[p2n, 2]; ξy += uy * Nij
            uz = grid.ps[p2n, 3]; ξz += uz * Nij
            # compute solid incremental deformation gradient
            df1 += ux * ∂Nx; df4 += uy * ∂Nx; df7 += uz * ∂Nx
            df2 += ux * ∂Ny; df5 += uy * ∂Ny; df8 += uz * ∂Ny
            df3 += ux * ∂Nz; df6 += uy * ∂Nz; df9 += uz * ∂Nz
        end)
        # update particle position
        mpts.ξ[ix, 1] += ξx
        mpts.ξ[ix, 2] += ξy
        mpts.ξ[ix, 3] += ξz
        # deformation gradient matrix
        update_F!(mpts, df1, df2, df3, df4, df5, df6, df7, df8, df9, ix)
        # update jacobian value and particle volume
        J = detF(mpts, ix)
        mpts.Ω[ix] = J * mpts.Ω0[ix]
        mpts.ρs[ix] = mpts.ρs0[ix] / J
        linearelastic!(mpts, df1, df2, df3, df4, df5, df6, df7, df8, df9, ix)
        druckerprager!(mpts, ix)
    end
end