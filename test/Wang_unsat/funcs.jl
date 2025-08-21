using KernelAbstractions
using Printf
using Symbolics

function tresetgridstatus!(grid::DeviceGrid{T1, T2}) where {T1, T2}
    fill!(grid.ms, T2(0.0))
    fill!(grid.fs, T2(0.0))
    fill!(grid.ps, T2(0.0))
    fill!(grid.vs, T2(0.0))
    fill!(grid.ext.mw, T2(0.0))
    fill!(grid.ext.mi, T2(0.0))
    fill!(grid.ext.fw, T2(0.0))
    fill!(grid.ext.pw, T2(0.0))
    fill!(grid.ext.vw, T2(0.0))
end

@kernel function tp2g!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        Ω, S, n, σw = mpts.Ω[ix], mpts.ext.S[ix], mpts.ext.n[ix], mpts.ext.σw[ix]
        ns, nl, σw = 1-n, S*n, mpts.ext.σw[ix]
        ms  = Ω * mpts.ρs[ix]
        mw  = Ω * mpts.ext.ρw
        mm  = ms * ns + mw * nl
        mmG = mm * mpts.G
        mwG = mw * mpts.G
        psx = ms * ns * mpts.vs[ix, 1]
        psy = ms * ns * mpts.vs[ix, 2]
        psz = ms * ns * mpts.vs[ix, 3]
        pwx = mw * nl * mpts.ext.vw[ix, 1]
        pwy = mw * nl * mpts.ext.vw[ix, 2]
        pwz = mw * nl * mpts.ext.vw[ix, 3]
        dragc = (nl * mw * 9.8) / (mpts.ext.kr[ix] * mpts.ext.k)
        dragx = dragc * (mpts.ext.vw[ix, 1] - mpts.vs[ix, 1])
        dragy = dragc * (mpts.ext.vw[ix, 2] - mpts.vs[ix, 2])
        dragz = dragc * (mpts.ext.vw[ix, 3] - mpts.vs[ix, 3])
        σxx, σyy, σzz = mpts.σij[ix, 1], mpts.σij[ix, 2], mpts.σij[ix, 3]
        σxy, σyz, σzx = mpts.σij[ix, 4], mpts.σij[ix, 5], mpts.σij[ix, 6]
        @KAunroll for iy in 1:mpts.NIC
            Nij = mpts.Nij[ix, iy]
            ∂Nx = mpts.∂Nx[ix, iy]
            ∂Ny = mpts.∂Ny[ix, iy]
            ∂Nz = mpts.∂Nz[ix, iy]
            p2n = mpts.p2n[ix, iy]
            # compute nodal mass
            @Σ grid.ext.mi[p2n] += Nij * mw
            @Σ grid.ext.mw[p2n] += Nij * mw * nl
            @Σ grid.ms[p2n] += Nij * ms * ns
            # compute nodal momentum
            @Σ grid.ext.vw[p2n, 1] += Nij * pwx
            @Σ grid.ext.vw[p2n, 2] += Nij * pwy
            @Σ grid.ext.vw[p2n, 3] += Nij * pwz
            @Σ grid.vs[p2n, 1] += Nij * psx
            @Σ grid.vs[p2n, 2] += Nij * psy
            @Σ grid.vs[p2n, 3] += Nij * psz
            # compute nodal total force
            @Σ grid.ext.fw[p2n, 1] += -∂Nx * Ω * σw - Nij * dragx
            @Σ grid.ext.fw[p2n, 2] += -∂Ny * Ω * σw - Nij * dragy
            @Σ grid.ext.fw[p2n, 3] += -∂Nz * Ω * σw - Nij * dragz + Nij * mwG
            @Σ grid.fs[p2n, 1] += -Ω * (∂Nx * (σxx + S*σw) + ∂Ny * σxy + ∂Nz * σzx)
            @Σ grid.fs[p2n, 2] += -Ω * (∂Ny * (σyy + S*σw) + ∂Nx * σxy + ∂Nz * σyz)
            @Σ grid.fs[p2n, 3] += -Ω * (∂Nz * (σzz + S*σw) + ∂Nx * σzx + ∂Ny * σyz) + Nij * mmG
        end
    end
end

@kernel function tsolvegrid!(grid::DeviceGrid{T1, T2}, Δt::T2) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        inv_mw = grid.ext.mw[ix] < eps(T2) ? T2(0.0) : 1 / grid.ext.mw[ix]
        inv_mi = grid.ext.mi[ix] < eps(T2) ? T2(0.0) : 1 / grid.ext.mi[ix]
        inv_ms = grid.ms[ix] < eps(T2) ? T2(0.0) : 1 / grid.ms[ix]
        # boundary condition
        grid.ext.vwxi[ix] ≠ T1(0) ? grid.ext.vw[ix, 1] = grid.ext.vwxv[ix] : nothing
        grid.ext.vwyi[ix] ≠ T1(0) ? grid.ext.vw[ix, 2] = grid.ext.vwyv[ix] : nothing
        grid.ext.vwzi[ix] ≠ T1(0) ? grid.ext.vw[ix, 3] = grid.ext.vwzv[ix] : nothing
        grid.vsxi[ix] ≠ T1(0) ? grid.vs[ix, 1] = grid.vsxv[ix] : nothing
        grid.vsyi[ix] ≠ T1(0) ? grid.vs[ix, 2] = grid.vsyv[ix] : nothing
        if grid.ext.topS[1] < 1.0
            grid.vszi[ix] ≠ T1(0) ? grid.vs[ix, 3] = grid.vszv[ix] : nothing
        end
        # compute nodal velocity
        grid.ext.vw[ix, 1] *= inv_mw
        grid.ext.vw[ix, 2] *= inv_mw
        grid.ext.vw[ix, 3] *= inv_mw
        grid.vs[ix, 1] *= inv_ms
        grid.vs[ix, 2] *= inv_ms
        grid.vs[ix, 3] *= inv_ms
        awx = grid.ext.fw[ix, 1] * inv_mi
        awy = grid.ext.fw[ix, 2] * inv_mi
        awz = grid.ext.fw[ix, 3] * inv_mi
        asx = (grid.fs[ix, 1] - grid.ext.mw[ix] * awx) * inv_ms
        asy = (grid.fs[ix, 2] - grid.ext.mw[ix] * awy) * inv_ms
        asz = (grid.fs[ix, 3] - grid.ext.mw[ix] * awz) * inv_ms
        # update nodal velocity
        grid.ext.vwT[ix, 1] = grid.ext.vw[ix, 1] + awx * Δt
        grid.ext.vwT[ix, 2] = grid.ext.vw[ix, 2] + awy * Δt
        grid.ext.vwT[ix, 3] = grid.ext.vw[ix, 3] + awz * Δt
        grid.vsT[ix, 1] = grid.vs[ix, 1] + asx * Δt
        grid.vsT[ix, 2] = grid.vs[ix, 2] + asy * Δt
        grid.vsT[ix, 3] = grid.vs[ix, 3] + asz * Δt
        # boundary condition
        grid.ext.vwxi[ix] ≠ T1(0) ? grid.ext.vwT[ix, 1] = grid.ext.vwxv[ix] : nothing
        grid.ext.vwyi[ix] ≠ T1(0) ? grid.ext.vwT[ix, 2] = grid.ext.vwyv[ix] : nothing
        if grid.ext.topS[1] < 1.0
            grid.ext.vwzi[ix] ≠ T1(0) ? grid.ext.vwT[ix, 3] = grid.ext.vwzv[ix] : nothing
        end
        grid.vsxi[ix] ≠ T1(0) ? grid.vsT[ix, 1] = grid.vsxv[ix] : nothing
        grid.vsyi[ix] ≠ T1(0) ? grid.vsT[ix, 2] = grid.vsyv[ix] : nothing
        grid.vszi[ix] ≠ T1(0) ? grid.vsT[ix, 3] = grid.vszv[ix] : nothing
    end
end

@kernel function tdoublemapping1!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}, Δt::T2) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        FLIP = mpts.FLIP
        ξsx = ξsy = ξsz = ξwx = ξwy = ξwz = T2(0.0)
        vsx = vsy = vsz = vwx = vwy = vwz = T2(0.0)
        @KAunroll for iy in 1:mpts.NIC
            Nij = mpts.Nij[ix, iy]
            p2n = mpts.p2n[ix, iy]
            ξwx += Nij * grid.ext.vwT[p2n, 1]
            ξwy += Nij * grid.ext.vwT[p2n, 2]
            ξwz += Nij * grid.ext.vwT[p2n, 3]
            ξsx += Nij * grid.vsT[p2n, 1]
            ξsy += Nij * grid.vsT[p2n, 2]
            ξsz += Nij * grid.vsT[p2n, 3]
            vwx += Nij * (grid.ext.vwT[p2n, 1] - grid.ext.vw[p2n, 1])
            vwy += Nij * (grid.ext.vwT[p2n, 2] - grid.ext.vw[p2n, 2])
            vwz += Nij * (grid.ext.vwT[p2n, 3] - grid.ext.vw[p2n, 3])
            vsx += Nij * (grid.vsT[p2n, 1] - grid.vs[p2n, 1])
            vsy += Nij * (grid.vsT[p2n, 2] - grid.vs[p2n, 2])
            vsz += Nij * (grid.vsT[p2n, 3] - grid.vs[p2n, 3])
        end
        # update particle position
        mpts.ξ[ix, 1] += Δt * ξsx
        mpts.ξ[ix, 2] += Δt * ξsy
        mpts.ξ[ix, 3] += Δt * ξsz
        # update particle velocity
        mpts.ext.vw[ix, 1] = FLIP * (mpts.ext.vw[ix, 1] + vwx) + (T2(1.0) - FLIP) * ξwx
        mpts.ext.vw[ix, 2] = FLIP * (mpts.ext.vw[ix, 2] + vwy) + (T2(1.0) - FLIP) * ξwy
        mpts.ext.vw[ix, 3] = FLIP * (mpts.ext.vw[ix, 3] + vwz) + (T2(1.0) - FLIP) * ξwz
        mpts.vs[ix, 1] = FLIP * (mpts.vs[ix, 1] + vsx) + (T2(1.0) - FLIP) * ξsx
        mpts.vs[ix, 2] = FLIP * (mpts.vs[ix, 2] + vsy) + (T2(1.0) - FLIP) * ξsy
        mpts.vs[ix, 3] = FLIP * (mpts.vs[ix, 3] + vsz) + (T2(1.0) - FLIP) * ξsz
    end
end

@kernel function tdoublemapping2!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        ms = mpts.ρs[ix] * mpts.Ω[ix] * (1.0 - mpts.ext.n[ix])
        mw = mpts.ext.ρw * mpts.Ω[ix] * mpts.ext.S[ix] * mpts.ext.n[ix]
        psx, psy, psz = ms * mpts.vs[ix, 1], ms * mpts.vs[ix, 2], ms * mpts.vs[ix, 3]
        pwx, pwy, pwz = mw * mpts.ext.vw[ix, 1], mw * mpts.ext.vw[ix, 2], mw * mpts.ext.vw[ix, 3]
        # update particle position & velocity
        @KAunroll for iy in 1:mpts.NIC
            Nij = mpts.Nij[ix, iy]
            p2n = mpts.p2n[ix, iy]
            @Σ grid.ext.pw[p2n, 1] += pwx * Nij
            @Σ grid.ext.pw[p2n, 2] += pwy * Nij
            @Σ grid.ext.pw[p2n, 3] += pwz * Nij
            @Σ grid.ps[p2n, 1] += psx * Nij
            @Σ grid.ps[p2n, 2] += psy * Nij
            @Σ grid.ps[p2n, 3] += psz * Nij
        end
    end
end

@kernel function tdoublemapping3!(grid::DeviceGrid{T1, T2}, Δt::T2) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        inv_mw = grid.ext.mw[ix] < eps(T2) ? T2(0.0) : 1 / grid.ext.mw[ix]
        inv_ms = grid.ms[ix] < eps(T2) ? T2(0.0) : 1 / grid.ms[ix]
        # compute nodal velocities
        grid.ext.pw[ix, 1] *= inv_mw * Δt
        grid.ext.pw[ix, 2] *= inv_mw * Δt
        grid.ext.pw[ix, 3] *= inv_mw * Δt
        grid.ps[ix, 1] *= inv_ms * Δt
        grid.ps[ix, 2] *= inv_ms * Δt
        grid.ps[ix, 3] *= inv_ms * Δt
        # fixed Dirichlet nodes
        grid.ext.vwxi[ix] ≠ T1(0) ? grid.ext.pw[ix, 1] = grid.ext.vwxv[ix] : nothing
        grid.ext.vwyi[ix] ≠ T1(0) ? grid.ext.pw[ix, 2] = grid.ext.vwyv[ix] : nothing
        if grid.ext.topS[1] < 1.0
            grid.ext.vwzi[ix] ≠ T1(0) ? grid.ext.pw[ix, 3] = grid.ext.vwzv[ix] : nothing
        end
        grid.vsxi[ix] ≠ T1(0) ? grid.ps[ix, 1] = grid.vsxv[ix] : nothing
        grid.vsyi[ix] ≠ T1(0) ? grid.ps[ix, 2] = grid.vsyv[ix] : nothing
        grid.vszi[ix] ≠ T1(0) ? grid.ps[ix, 3] = grid.vszv[ix] : nothing
    end
end

# let
#     @variables x
#     y = 1-0.10152*(x/9800)^2.4279
#     dy_dx = Symbolics.derivative(y, x; simplify=true)
#     println(dy_dx)
# end

∂S∂ψm(ψm) = -5.0290345007256016e-11 * (ψm^1.4279000000000002)
Sr(ψm) = 1 - 0.10152 * (ψm / 9800)^2.4279
kr(S) = 1 - 2.207 * (1 - S)

@kernel function tg2p!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}, material::Material, t_eld::T2, t_cur::T2, Δt::T2) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mpts.np
        dfs1 = dfs2 = dfs3 = dfs4 = dfs5 = dfs6 = dfs7 = dfs8 = dfs9 = T2(0.0)
        dfw1 = dfw2 = dfw3 = dfw4 = dfw5 = dfw6 = dfw7 = dfw8 = dfw9 = T2(0.0)
        @KAunroll for iy in 1:mpts.NIC
            p2n = mpts.p2n[ix, iy]
            ∂Nx = mpts.∂Nx[ix, iy]; vsx = grid.ps[p2n, 1]; vwx = grid.ext.pw[p2n, 1]
            ∂Ny = mpts.∂Ny[ix, iy]; vsy = grid.ps[p2n, 2]; vwy = grid.ext.pw[p2n, 2]
            ∂Nz = mpts.∂Nz[ix, iy]; vsz = grid.ps[p2n, 3]; vwz = grid.ext.pw[p2n, 3]
            # compute solid incremental deformation gradient
            dfs1 += vsx * ∂Nx; dfs4 += vsy * ∂Nx; dfs7 += vsz * ∂Nx
            dfs2 += vsx * ∂Ny; dfs5 += vsy * ∂Ny; dfs8 += vsz * ∂Ny
            dfs3 += vsx * ∂Nz; dfs6 += vsy * ∂Nz; dfs9 += vsz * ∂Nz
            dfw1 += vwx * ∂Nx; dfw4 += vwy * ∂Nx; dfw7 += vwz * ∂Nx
            dfw2 += vwx * ∂Ny; dfw5 += vwy * ∂Ny; dfw8 += vwz * ∂Ny
            dfw3 += vwx * ∂Nz; dfw6 += vwy * ∂Nz; dfw9 += vwz * ∂Nz
        end
        # deformation gradient matrix
        update_F!(mpts, dfs1, dfs2, dfs3, dfs4, dfs5, dfs6, dfs7, dfs8, dfs9, ix)
        # update jacobian value and particle volume
        J = detF(mpts, ix)
        ΔJ = J * mpts.Ω0[ix] / mpts.Ω[ix]
        mpts.Ω[ix] = J * mpts.Ω0[ix]
        mpts.ρs[ix] = mpts.ρs[ix] / ΔJ
        # update pore pressure
        ψm, S, n = -mpts.ext.σw[ix], mpts.ext.S[ix], mpts.ext.n[ix]
        if ψm > T2(0.0)
            #@info "ψm = $(ψm), S = $(S), n = $(n)"
            Cs = ∂S∂ψm(ψm)
            # update hydraulic conductivity and the degree of saturation
            mpts.ext.S[ix] = Sr(-mpts.ext.σw[ix])
            mpts.ext.kr[ix] = kr(mpts.ext.S[ix])
        else
            Cs = T2(0.0)
            S = T2(1.0)
            # update hydraulic conductivity and the degree of saturation
            mpts.ext.S[ix] = T2(1.0)
            mpts.ext.kr[ix] = T2(1.0)
        end
        mpts.ext.σw[ix] += (1 / (n*Cs + (n*S)/mpts.ext.Kw)) * (
            (1.0 - n) * S * (dfs1 + dfs5 + dfs9) + 
                   n  * S * (dfw1 + dfw5 + dfw9))
        # update porosity
        mpts.ext.n[ix] = clamp(1.0 - (1.0 - mpts.ext.n[ix]) / ΔJ, 0.0, 1.0)
        material!(mpts, t_eld, t_cur, Δt, 
            dfs1, dfs2, dfs3, dfs4, dfs5, dfs6, dfs7, dfs8, dfs9, 
            ix, material)
    end
end

function tprocedure!(conf::Config, grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}, vpos, vcol) where {T1, T2}
    model_info(conf, grid, mpts)
    t_cur = T2(conf.t_cur)
    t_tol = T2(conf.t_tol)
    t_eld = T2(conf.t_eld)
    Δt    = T2(conf.Δt)
    dev   = conf.dev
    h5    = conf.h5
    dev_grid, dev_mpts = host2device(dev, grid, mpts)

    fid = set_hdf5(conf)
    printer = set_pb(conf)
    while t_cur < t_tol
        hdf5!(h5, fid, t_cur, mpts, dev_mpts)

        vpos[] = dev_mpts.ξ
        vcol[] = dev_mpts.ext.vw[:, 3]
        tresetgridstatus!(dev_grid)
        resetmpstatus!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, conf.basis)
        tp2g!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts)
        tsolvegrid!(dev)(ndrange=dev_grid.ni, dev_grid, Δt)
        tdoublemapping1!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, Δt)
        tdoublemapping2!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts)
        tdoublemapping3!(dev)(ndrange=dev_grid.ni, dev_grid, Δt)
        tg2p!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, conf.material, t_eld, t_cur, Δt)
        dev_mpts.ext.σw[1:8] .= 0.0
        dev_mpts.ext.S[1:8] .= 1.0
        
        mpts_top_idx = findall(i->dev_mpts.ξ[i, 3]≥1-grid.h, 1:mpts.np)
        mpts_top_inp = length(mpts_top_idx)
        dev_grid.ext.topS[1] = sum(dev_mpts.ext.S[mpts_top_idx]) / mpts_top_inp

        t_cur += Δt
        h5.iters[] += 1
        update_pb!(printer, t_cur, t_tol)

        sleep(0.5)
    end
    finish_pb!(conf, printer); KAsync(dev)
    device2host!(mpts, dev_mpts)
    hdf5!(h5, fid, grid)
    close(fid)
end