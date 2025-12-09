function resetgridstatus!(grid::DeviceGrid{T1, T2}) where {T1, T2}
    fill!(grid.ext.ms, T2(0.0))
    fill!(grid.ext.fs, T2(0.0))
    fill!(grid.ext.ps, T2(0.0))
    fill!(grid.ext.vs, T2(0.0))
end

@kernel function p2g!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}, Gg::T2) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        Ω, ms = mpts.ext.Ω[ix], mpts.ext.ρs[ix] * mpts.ext.Ω[ix]
        msG = ms * Gg
        psx, psy, psz = mpts.ext.vs[ix, 1] * ms, mpts.ext.vs[ix, 2] * ms, mpts.ext.vs[ix, 3] * ms
        σxx, σyy, σzz = mpts.ext.σij[ix, 1], mpts.ext.σij[ix, 2], mpts.ext.σij[ix, 3]
        σxy, σyz, σzx = mpts.ext.σij[ix, 4], mpts.ext.σij[ix, 5], mpts.ext.σij[ix, 6]
        @N∂NBspline2(grid, mpts, ix, begin
            # compute nodal mass
            @Σ grid.ext.ms[p2n] += Nij * ms
            # compute nodal momentum
            @Σ grid.ext.vs[p2n, 1] += Nij * psx
            @Σ grid.ext.vs[p2n, 2] += Nij * psy
            @Σ grid.ext.vs[p2n, 3] += Nij * psz
            # compute nodal total force for solid
            @Σ grid.ext.fs[p2n, 1] += -Ω * (∂Nx * σxx + ∂Ny * σxy + ∂Nz * σzx)
            @Σ grid.ext.fs[p2n, 2] += -Ω * (∂Nx * σxy + ∂Ny * σyy + ∂Nz * σyz)
            @Σ grid.ext.fs[p2n, 3] += -Ω * (∂Nx * σzx + ∂Ny * σyz + ∂Nz * σzz) + Nij * msG
        end)
    end
end

@kernel function solvegrid!(grid::DeviceGrid{T1, T2}, Δt::T2) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni && grid.ext.ms[ix] ≠ T2(0.0)
        ms_denom = grid.ext.ms[ix] < eps(T2) ? T2(0.0) : 1 / grid.ext.ms[ix]
        # boundary condition
        gx, gy, gz = get_grid_ξ(grid, ix)
        bcx = gx ≤ T2(0.0) || gx ≥ T2(0.05) || gz ≤ T2(0.0) || gy ≤ T2(0.0)
        bcy = gy ≤ T2(0.0) || gz ≤ T2(0.00)
        bcz = gz ≤ T2(0.0)
        grid.ext.vs[ix, 1] = bcx ? T2(0.0) : grid.ext.vs[ix, 1]
        grid.ext.vs[ix, 2] = bcy ? T2(0.0) : grid.ext.vs[ix, 2]
        grid.ext.vs[ix, 3] = bcz ? T2(0.0) : grid.ext.vs[ix, 3]
        # compute nodal velocity
        grid.ext.vs[ix, 1] *= ms_denom
        grid.ext.vs[ix, 2] *= ms_denom
        grid.ext.vs[ix, 3] *= ms_denom
        # damping force for solid
        dampvs = -grid.ζs * sqrt(grid.ext.fs[ix, 1] * grid.ext.fs[ix, 1] + 
                                 grid.ext.fs[ix, 2] * grid.ext.fs[ix, 2] + 
                                 grid.ext.fs[ix, 3] * grid.ext.fs[ix, 3])
        # compute nodal total force for mixture
        Fs_x = grid.ext.fs[ix, 1] + dampvs * sign(grid.ext.vs[ix, 1])
        Fs_y = grid.ext.fs[ix, 2] + dampvs * sign(grid.ext.vs[ix, 2])
        Fs_z = grid.ext.fs[ix, 3] + dampvs * sign(grid.ext.vs[ix, 3])
        # update nodal velocity
        grid.ext.vsT[ix, 1] = grid.ext.vs[ix, 1] + Fs_x * Δt * ms_denom
        grid.ext.vsT[ix, 2] = grid.ext.vs[ix, 2] + Fs_y * Δt * ms_denom
        grid.ext.vsT[ix, 3] = grid.ext.vs[ix, 3] + Fs_z * Δt * ms_denom
        # boundary condition
        grid.ext.vsT[ix, 1] = bcx ? T2(0.0) : grid.ext.vsT[ix, 1]
        grid.ext.vsT[ix, 2] = bcy ? T2(0.0) : grid.ext.vsT[ix, 2]
        grid.ext.vsT[ix, 3] = bcz ? T2(0.0) : grid.ext.vsT[ix, 3]
        # compute velocity increment
        grid.ext.vs[ix, 1] = grid.ext.vsT[ix, 1] - grid.ext.vs[ix, 1]
        grid.ext.vs[ix, 2] = grid.ext.vsT[ix, 2] - grid.ext.vs[ix, 2]
        grid.ext.vs[ix, 3] = grid.ext.vsT[ix, 3] - grid.ext.vs[ix, 3]
    end
end

@kernel function doublemapping1!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        FLIP = mpts.FLIP
        ξx = ξy = ξz = vx = vy = vz = T2(0.0)
        @NijBspline2(grid, mpts, ix, begin
            ξx += Nij * grid.ext.vsT[p2n, 1]
            ξy += Nij * grid.ext.vsT[p2n, 2]
            ξz += Nij * grid.ext.vsT[p2n, 3]
            vx += Nij * grid.ext.vs[p2n, 1]
            vy += Nij * grid.ext.vs[p2n, 2]
            vz += Nij * grid.ext.vs[p2n, 3]
        end)
        # update particle velocity
        mpts.ext.vs[ix, 1] = FLIP * (mpts.ext.vs[ix, 1] + vx) + (T2(1.0) - FLIP) * ξx
        mpts.ext.vs[ix, 2] = FLIP * (mpts.ext.vs[ix, 2] + vy) + (T2(1.0) - FLIP) * ξy
        mpts.ext.vs[ix, 3] = FLIP * (mpts.ext.vs[ix, 3] + vz) + (T2(1.0) - FLIP) * ξz
        # update grid momentum
        ms = mpts.ext.ρs[ix] * mpts.ext.Ω[ix]
        px, py, pz = ms * mpts.ext.vs[ix, 1], ms * mpts.ext.vs[ix, 2], ms * mpts.ext.vs[ix, 3]
        # update CFL conditions
        nid  = mpts.nid[ix]
        Ks   = mpts.ext.Ks[nid]
        Gs   = mpts.ext.Gs[nid]
        cdil = sqrt((Ks + T2(1.333333) * Gs) / mpts.ext.ρs[ix]) # 4/3 ≈ 1.333333
        mpts.cfl[ix] = min(grid.h / (cdil + abs(mpts.ext.vs[ix, 1])), 
                           grid.h / (cdil + abs(mpts.ext.vs[ix, 2]))) 
        # update grid momentum
        @NijBspline2(grid, mpts, ix, begin
            @Σ grid.ext.ps[p2n, 1] += px * Nij
            @Σ grid.ext.ps[p2n, 2] += py * Nij
            @Σ grid.ext.ps[p2n, 3] += pz * Nij
        end)
    end
end

@kernel function doublemapping2!(grid::DeviceGrid{T1, T2}, Δt::T2) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni && grid.ext.ms[ix] ≠ T2(0.0)
        ms_denom = grid.ext.ms[ix] < eps(T2) ? T2(0.0) : 1 / grid.ext.ms[ix]
        # compute nodal velocities
        grid.ext.ps[ix, 1] = grid.ext.ps[ix, 1] * ms_denom * Δt
        grid.ext.ps[ix, 2] = grid.ext.ps[ix, 2] * ms_denom * Δt
        grid.ext.ps[ix, 3] = grid.ext.ps[ix, 3] * ms_denom * Δt
        # fixed Dirichlet nodes
        gx, gy, gz = get_grid_ξ(grid, ix)
        bcx = gx ≤ T2(0.0) || gx ≥ T2(0.05) || gz ≤ T2(0.0) || gy ≤ T2(0.0)
        bcy = gy ≤ T2(0.0) || gz ≤ T2(0.00)
        bcz = gz ≤ T2(0.0)
        grid.ext.ps[ix, 1] = bcx ? T2(0.0) : grid.ext.ps[ix, 1]
        grid.ext.ps[ix, 2] = bcy ? T2(0.0) : grid.ext.ps[ix, 2]
        grid.ext.ps[ix, 3] = bcz ? T2(0.0) : grid.ext.ps[ix, 3]
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
            ux = grid.ext.ps[p2n, 1]; ξx += ux * Nij
            uy = grid.ext.ps[p2n, 2]; ξy += uy * Nij
            uz = grid.ext.ps[p2n, 3]; ξz += uz * Nij
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
        mpts.ext.Ω[ix] = J * mpts.ext.Ω0
        mpts.ext.ρs[ix] = mpts.ext.ρs0 / J
        linearelastic!(mpts, df1, df2, df3, df4, df5, df6, df7, df8, df9, ix)
        druckerprager!(mpts, ix)
    end
end