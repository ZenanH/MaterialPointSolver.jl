using KernelAbstractions

@kernel function tp2g!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}, G::T2) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        Ω, ms = mpts.Ω[ix], mpts.ρs[ix] * mpts.Ω[ix]
        msG = ms * G
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

@kernel function tsolvegrid!(grid::DeviceGrid{T1, T2}, Δt::T2) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : 1 / grid.ms[ix]
        # boundary condition
        gx, gy, gz = _get_grid_coord(grid, ix)
        bc = gx ≤ grid.x1 + grid.h * T2(6.0) || gx ≥ grid.x2 - grid.h * T2(6.0) ||
             gy ≤ grid.y1 + grid.h * T2(6.0) || gy ≥ grid.y2 - grid.h * T2(6.0) ||
             gz ≤ grid.z1 + grid.h * T2(6.0)# || gz ≥ grid.z2 - grid.h * T2(6.0)
        if bc
            grid.vs[ix, 1] = T2(0.0)
            grid.vs[ix, 2] = T2(0.0)
            grid.vs[ix, 3] = T2(0.0)
        end
        # compute nodal velocity
        grid.vs[ix, 1] *= ms_denom
        grid.vs[ix, 2] *= ms_denom
        grid.vs[ix, 3] *= ms_denom
        # update nodal force
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
        # add friction boundary condition
        nv_id = grid.ext.id[ix]
        if nv_id ≠ T1(0)
            μ = grid.ext.μ
            nx, ny, nz = grid.ext.n[nv_id, 1], grid.ext.n[nv_id, 2], grid.ext.n[nv_id, 3]
            vx, vy, vz = grid.vsT[ix, 1], grid.vsT[ix, 2], grid.vsT[ix, 3]
            m = grid.ms[ix]
            v_n = vx * nx + vy * ny + vz * nz
            if v_n < T2(0.0)
                J_n     = -m * v_n
                v_tx    = vx - v_n * nx
                v_ty    = vy - v_n * ny
                v_tz    = vz - v_n * nz
                v_t     = sqrt(v_tx * v_tx + v_ty * v_ty + v_tz * v_tz)
                J_t_max = μ * abs(J_n)
                J_tx, J_ty, J_tz = T2(0.0), T2(0.0), T2(0.0)
                if v_t > T2(0.0)
                    J_t_required = m * v_t    
                    if J_t_required > J_t_max
                        scale = J_t_max / v_t
                        J_tx  = -scale * v_tx
                        J_ty  = -scale * v_ty
                        J_tz  = -scale * v_tz
                    else
                        J_tx = -m * v_tx
                        J_ty = -m * v_ty
                        J_tz = -m * v_tz
                    end                    
                end
                grid.vsT[ix, 1] = vx + (J_n * nx + J_tx) * ms_denom
                grid.vsT[ix, 2] = vy + (J_n * ny + J_ty) * ms_denom
                grid.vsT[ix, 3] = vz + (J_n * nz + J_tz) * ms_denom
            end
        end
        # boundary condition
        if bc
            grid.vsT[ix, 1] = T2(0.0)
            grid.vsT[ix, 2] = T2(0.0)
            grid.vsT[ix, 3] = T2(0.0)
        end
    end
end

@kernel function tdoublemapping1!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}, Δt::T2) where {T1, T2}
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
        # update tmp particle velocity
        vx = FLIP * (mpts.vs[ix, 1] + vx) + (T2(1.0) - FLIP) * ξx
        vy = FLIP * (mpts.vs[ix, 2] + vy) + (T2(1.0) - FLIP) * ξy
        vz = FLIP * (mpts.vs[ix, 3] + vz) + (T2(1.0) - FLIP) * ξz
        # update tmp particle position
        x = mpts.ξ[ix, 1] + Δt * ξx
        y = mpts.ξ[ix, 2] + Δt * ξy
        z = mpts.ξ[ix, 3] + Δt * ξz
        # find the nearest node
        bid, x0, y0, z0 = _get_grid_coord(grid, x, y, z)
        # correct particle position and velocity (friction boundary)
        nv_id = grid.ext.id[bid]
        if nv_id ≠ T1(0)
            nx, ny, nz = grid.ext.n[nv_id, 1], grid.ext.n[nv_id, 2], grid.ext.n[nv_id, 3]
            d  = (x - x0) * nx + (y - y0) * ny + (z - z0) * nz
            if d < T2(0.0)
                x -= d * nx
                y -= d * ny
                z -= d * nz
                v_n = vx * nx + vy * ny + vz * nz
                if v_n < T2(0.0)
                    vx -= v_n * nx
                    vy -= v_n * ny
                    vz -= v_n * nz
                end
            end
        end
        # update particle position
        mpts.ξ[ix, 1] = x
        mpts.ξ[ix, 2] = y
        mpts.ξ[ix, 3] = z
        # update particle velocity
        mpts.vs[ix, 1] = vx
        mpts.vs[ix, 2] = vy
        mpts.vs[ix, 3] = vz
    end
end

@kernel function tdoublemapping3!(grid::DeviceGrid{T1, T2}, Δt::T2) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : 1 / grid.ms[ix]
        # compute nodal velocities
        grid.ps[ix, 1] = grid.ps[ix, 1] * ms_denom * Δt
        grid.ps[ix, 2] = grid.ps[ix, 2] * ms_denom * Δt
        grid.ps[ix, 3] = grid.ps[ix, 3] * ms_denom * Δt
        # fixed Dirichlet nodes
        gx, gy, gz = _get_grid_coord(grid, ix)
        if gx ≤ grid.x1 + grid.h * T2(6.0) || gx ≥ grid.x2 - grid.h * T2(6.0) ||
           gy ≤ grid.y1 + grid.h * T2(6.0) || gy ≥ grid.y2 - grid.h * T2(6.0) ||
           gz ≤ grid.z1 + grid.h * T2(6.0)# || gz ≥ grid.z2 - grid.h * T2(6.0)
            grid.ps[ix, 1] = T2(0.0)
            grid.ps[ix, 2] = T2(0.0)
            grid.ps[ix, 3] = T2(0.0)
        end
    end
end

@kernel function tg2p!(
    grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}, Δt::T2
) where {T1, T2}
    ix = @index(Global)
    invΔt = T2(1.0) / Δt
    if ix ≤ mpts.np
        ϵ0̇ = mpts.ext.ϵ0̇
        A  = mpts.ext.A
        invn = T2(0.333333)
        Ks   = mpts.Ks[mpts.nid[ix]]
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

        D1 = df1 * invΔt
        D5 = df5 * invΔt
        D9 = df9 * invΔt
        D2 = T2(0.5) * (df2 + df4) * invΔt
        D3 = T2(0.5) * (df3 + df7) * invΔt
        D6 = T2(0.5) * (df6 + df8) * invΔt

        trD3 = T2(0.333333) * (D1 + D5 + D9)
        D1 -= trD3; D5 -= trD3; D9 -= trD3

        eqϵ̇ = sqrt(ϵ0̇^2 + T2(0.5) * ((D1^2 + D5^2 + D9^2) + T2(2.0) * (D2^2 + D3^2 + D6^2))) # equivalent strain rate

        η = T2(0.5) * A^(-invn) * eqϵ̇^(invn - 1) # effective viscosity

        τxx = T2(2.0) * η * D1
        τyy = T2(2.0) * η * D5
        τzz = T2(2.0) * η * D9
        τxy = T2(2.0) * η * D2
        τyz = T2(2.0) * η * D6
        τzx = T2(2.0) * η * D3

        P = T2(-0.333333) * (mpts.σij[ix, 1] + mpts.σij[ix, 2] + mpts.σij[ix, 3]) # mean stress

        P -= Ks * Δt * trD3 * T2(3.0) # pressure update

        mpts.σij[ix, 1] = -P + τxx
        mpts.σij[ix, 2] = -P + τyy
        mpts.σij[ix, 3] = -P + τzz
        mpts.σij[ix, 4] = τxy
        mpts.σij[ix, 5] = τyz
        mpts.σij[ix, 6] = τzx

        # deformation gradient matrix
        update_F!(mpts, df1, df2, df3, df4, df5, df6, df7, df8, df9, ix)
        # update jacobian value and particle volume
        J = detF(mpts, ix)
        mpts.Ω[ix] = J * mpts.Ω0[ix]
        mpts.ρs[ix] = mpts.ρs0[ix] / J
        mpts.ext.σm[ix] = -P
    end
end

function tprocedure!(conf::Config, grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
    model_info(conf, grid, mpts)
    t_cur = T2(conf.t_cur)
    t_tol = T2(conf.t_tol)
    t_eld = T2(conf.t_eld)
    Δt    = T2(conf.Δt)
    Gg    = T2(mpts.G)
    dev   = conf.dev
    h5    = conf.h5
    dev_grid, dev_mpts = host2device(dev, grid, mpts)

    fid = set_hdf5(conf)
    printer = set_pb(conf)
    while t_cur < t_tol
        G = t_cur ≤ t_eld ? (Gg * t_cur) / t_eld : Gg
        hdf5!(h5, fid, t_cur, mpts, dev_mpts)
        resetgridstatus!(dev_grid)
        resetmpstatus!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, conf.basis)
        tp2g!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, G)
        tsolvegrid!(dev)(ndrange=dev_grid.ni, dev_grid, Δt)
        tdoublemapping1!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, Δt)
        doublemapping2!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts)
        tdoublemapping3!(dev)(ndrange=dev_grid.ni, dev_grid, Δt)
        tg2p!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, Δt)
        t_cur += Δt
        h5.iters[] += 1
        update_pb!(printer, t_cur, t_tol)
    end
    finish_pb!(conf, printer); KAsync(dev)
    device2host!(mpts, dev_mpts)
    hdf5!(h5, fid, grid)
    close(fid)
end