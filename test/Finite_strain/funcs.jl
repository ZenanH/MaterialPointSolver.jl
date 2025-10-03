using KernelAbstractions
include(joinpath(@__DIR__, "utils.jl"))

function tresetgridstatus!(grid::DeviceGrid{T1, T2}) where {T1, T2}
    fill!(grid.ms, T2(0.0))
    fill!(grid.fs, T2(0.0))
    fill!(grid.ps, T2(0.0))
    fill!(grid.vs, T2(0.0))
end

@kernel function tresetmpstatus!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}, ::Bspline2Basis) where {T1, T2}
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

@kernel function tp2g!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
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

@kernel function tsolvegrid!(grid::DeviceGrid{T1, T2}, Δt::T2) where {T1, T2}
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

@kernel function tdoublemapping2!(grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
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

@kernel function tdoublemapping3!(grid::DeviceGrid{T1, T2}, Δt::T2) where {T1, T2}
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

@kernel function tg2p!(
    grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}, 
    material::Material, t_eld::T2, t_cur::T2, Δt::T2
) where {T1, T2}
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
        Jc = J > eps(T2) ? J : eps(T2)
        mpts.Ω[ix] = Jc * mpts.Ω0[ix]
        mpts.ρs[ix] = mpts.ρs0[ix] / Jc
        tliKE!(mpts, Jc, ix)
        #tdpP!(mpts, Jc, ix)
        tupdateHP!(mpts, Jc, ix)
    end
end

# Hencky 试算（使用 h_e = h - h_p）→ τ → σ，零分配/GPU OK
@inline Base.@propagate_inbounds function tliKE!(
    mpts::DeviceParticle{T1,T2}, J::T2, ix::Int
) where {T1,T2}
    nid = mpts.nid[ix]
    Gs  = mpts.Gs[nid]   # μ
    Ks  = mpts.Ks[nid]   # κ

    # F
    F11 = mpts.F[ix, 1]; F12 = mpts.F[ix, 2]; F13 = mpts.F[ix, 3]
    F21 = mpts.F[ix, 4]; F22 = mpts.F[ix, 5]; F23 = mpts.F[ix, 6]
    F31 = mpts.F[ix, 7]; F32 = mpts.F[ix, 8]; F33 = mpts.F[ix, 9]

    # 夹紧后的 J
    J = J > eps(T2) ? J : eps(T2)
    invJ = one(T2) / J

    # C = FᵀF
    C11 = F11*F11 + F21*F21 + F31*F31
    C12 = F11*F12 + F21*F22 + F31*F32
    C13 = F11*F13 + F21*F23 + F31*F33
    C22 = F12*F12 + F22*F22 + F32*F32
    C23 = F12*F13 + F22*F23 + F32*F33
    C33 = F13*F13 + F23*F23 + F33*F33

    b1,b2,b3, q11,q12,q13, q21,q22,q23, q31,q32,q33 =
        sym_eig3x3_jacobi(C11, C22, C33, C12, C13, C23, 6)

    λ1 = b1 > T2(1e-12) ? b1 : T2(1e-12)
    λ2 = b2 > T2(1e-12) ? b2 : T2(1e-12)
    λ3 = b3 > T2(1e-12) ? b3 : T2(1e-12)
    d1 = T2(0.5) * log(λ1)
    d2 = T2(0.5) * log(λ2)
    d3 = T2(0.5) * log(λ3)

    # 总对数应变 h（并写回保存，供 DP 后回推 hp 用）
    h11 = d1*q11*q11 + d2*q12*q12 + d3*q13*q13
    h22 = d1*q21*q21 + d2*q22*q22 + d3*q23*q23
    h33 = d1*q31*q31 + d2*q32*q32 + d3*q33*q33
    h12 = d1*q11*q21 + d2*q12*q22 + d3*q13*q23
    h13 = d1*q11*q31 + d2*q12*q32 + d3*q13*q33
    h23 = d1*q21*q31 + d2*q22*q32 + d3*q23*q33

    mpts.ext.hij[ix, 1] = h11; mpts.ext.hij[ix, 2] = h22; mpts.ext.hij[ix, 3] = h33
    mpts.ext.hij[ix, 4] = h23; mpts.ext.hij[ix, 5] = h13; mpts.ext.hij[ix, 6] = h12

    # 取历史塑性对数应变 h_p
    hp11 = mpts.ext.hpij[ix, 1]; hp22 = mpts.ext.hpij[ix, 2]; hp33 = mpts.ext.hpij[ix, 3]
    hp23 = mpts.ext.hpij[ix, 4]; hp13 = mpts.ext.hpij[ix, 5]; hp12 = mpts.ext.hpij[ix, 6]

    # 弹性对数应变 h_e = h - h_p
    he11 = h11 - hp11; he22 = h22 - hp22; he33 = h33 - hp33
    he12 = h12 - hp12; he13 = h13 - hp13; he23 = h23 - hp23

    trh = he11 + he22 + he33
    dev11 = he11 - T2(0.333333)*trh
    dev22 = he22 - T2(0.333333)*trh
    dev33 = he33 - T2(0.333333)*trh
    dev12 = he12; dev13 = he13; dev23 = he23

    # Kirchhoff 试算
    τ11 = Ks*trh + (T2(2)*Gs)*dev11
    τ22 = Ks*trh + (T2(2)*Gs)*dev22
    τ33 = Ks*trh + (T2(2)*Gs)*dev33
    τ12 = (T2(2)*Gs)*dev12
    τ13 = (T2(2)*Gs)*dev13
    τ23 = (T2(2)*Gs)*dev23

    # Cauchy
    mpts.σij[ix,1] = τ11 * invJ
    mpts.σij[ix,2] = τ22 * invJ
    mpts.σij[ix,3] = τ33 * invJ
    mpts.σij[ix,4] = τ12 * invJ
    mpts.σij[ix,5] = τ23 * invJ
    mpts.σij[ix,6] = τ13 * invJ
end


@inline Base.@propagate_inbounds function tdpP!(mpts::DeviceParticle{T1, T2}, J::T2, ix::Int) where {T1, T2}
    nid = mpts.nid[ix]
    c   = mpts.c[nid]
    ϕ   = mpts.ϕ[nid]
    ψ   = mpts.ψ[nid]
    σt  = mpts.σt[nid]
    Gs  = mpts.Gs[nid]
    Ks  = mpts.Ks[nid]

    # update deviatoric stress tensor
    σm = (mpts.σij[ix, 1] + mpts.σij[ix, 2] + mpts.σij[ix, 3]) * T2(0.333333)
    s1 = mpts.σij[ix, 1] - σm
    s2 = mpts.σij[ix, 2] - σm
    s3 = mpts.σij[ix, 3] - σm
    s4 = mpts.σij[ix, 4]
    s5 = mpts.σij[ix, 5]
    s6 = mpts.σij[ix, 6]

    # drucker-prager
    τ  = sqrt(T2(0.5) * (s1 * s1 + s2 * s2 + s3 * s3) +
                         s4 * s4 + s5 * s5 + s6 * s6)
    τ  = τ < eps(T2) ? eps(T2) : τ
    kϕ = (T2(6.0) * c * cos(ϕ)) / (T2(1.732051) * (T2(3.0) + sin(ϕ)))
    qϕ = (T2(6.0) *     sin(ϕ)) / (T2(1.732051) * (T2(3.0) + sin(ϕ)))
    qψ = (T2(6.0) *     sin(ψ)) / (T2(1.732051) * (T2(3.0) + sin(ψ)))
    σt = min(σt, kϕ / qϕ)
    αb = sqrt(T2(1.0) + qϕ * qϕ) - qϕ
    τb = kϕ - qϕ * σt
    fs = τ + qϕ * σm - kϕ            # yield function considering shear failure
    ft = σm - σt                     # yield function considering tensile failure
    BF = (τ - τb) - (αb * (σm - σt)) # BF is used to classify shear failure from tensile failure

    # determination of failure criteria
    ## shear failure correction
    if ((σm < σt) && (fs > T2(0.0))) ||
       ((σm ≥ σt) && (BF > T2(0.0)))
        Δλs  = (J * fs) / (Gs + Ks * qϕ * qψ) # *** 修改 ***：Δλs 乘以 J（σ 空间等效刚度为原来的 1/J）
        tmp1 = σm - (Ks / J) * qψ * Δλs # *** 修改 ***：体积项 Ks → Ks/J
        tmp2 = (kϕ - qϕ * tmp1) / τ
        mpts.σij[ix, 1] = s1 * tmp2 + tmp1
        mpts.σij[ix, 2] = s2 * tmp2 + tmp1
        mpts.σij[ix, 3] = s3 * tmp2 + tmp1
        mpts.σij[ix, 4] = s4 * tmp2
        mpts.σij[ix, 5] = s5 * tmp2
        mpts.σij[ix, 6] = s6 * tmp2
        mpts.ϵq[ix] += (Δλs / J) * sqrt(T2(0.333333) + T2(0.222222) * qψ * qψ)
        mpts.ϵk[ix] += (Δλs / J) * qψ
    end

    ## tensile failure correction
    if (σm ≥ σt) && (BF ≤ T2(0.0))
        Δλt = (J * ft) / Ks # *** 修改 ***：Δλt 乘以 J
        mpts.σij[ix, 1] = s1 + σt
        mpts.σij[ix, 2] = s2 + σt
        mpts.σij[ix, 3] = s3 + σt
        mpts.ϵq[ix] += (Δλt / J) * T2(0.333333) * T2(1.414214)
        mpts.ϵk[ix] += (Δλt / J)
    end
end

# 在 DP 更新 σ 后：τ = J*σ → 反解 h_e → h_p = h - h_e，零分配
@inline Base.@propagate_inbounds function tupdateHP!(
    mpts::DeviceParticle{T1,T2}, J::T2, ix::Int
) where {T1, T2}
    nid = mpts.nid[ix]
    Gs  = mpts.Gs[nid]   # μ
    Ks  = mpts.Ks[nid]   # κ

    # τ = J σ
    τ11 = J * mpts.σij[ix,1]
    τ22 = J * mpts.σij[ix,2]
    τ33 = J * mpts.σij[ix,3]
    τ23 = J * mpts.σij[ix,4]
    τ13 = J * mpts.σij[ix,5]
    τ12 = J * mpts.σij[ix,6]

    trτ = τ11 + τ22 + τ33

    # dev τ
    s11 = τ11 - T2(0.333333)*trτ
    s22 = τ22 - T2(0.333333)*trτ
    s33 = τ33 - T2(0.333333)*trτ
    s23 = τ23
    s13 = τ13
    s12 = τ12

    # 反解 h_e = (1/(3K)) tr(τ) I + (1/(2μ)) dev(τ)
    inv3K = one(T2) / (T2(3)*Ks)
    inv2G = one(T2) / (T2(2)*Gs)
    he11 = trτ*inv3K + s11*inv2G
    he22 = trτ*inv3K + s22*inv2G
    he33 = trτ*inv3K + s33*inv2G
    he23 = s23*inv2G
    he13 = s13*inv2G
    he12 = s12*inv2G

    # h_p = h - h_e（h 已在 tliKE_with_hp! 里存到 hij）
    h11 = mpts.ext.hij[ix, 1]; h22 = mpts.ext.hij[ix, 2]; h33 = mpts.ext.hij[ix, 3]
    h23 = mpts.ext.hij[ix, 4]; h13 = mpts.ext.hij[ix, 5]; h12 = mpts.ext.hij[ix, 6]

    mpts.ext.hpij[ix, 1] = h11 - he11
    mpts.ext.hpij[ix, 2] = h22 - he22
    mpts.ext.hpij[ix, 3] = h33 - he33
    mpts.ext.hpij[ix, 4] = h23 - he23
    mpts.ext.hpij[ix, 5] = h13 - he13
    mpts.ext.hpij[ix, 6] = h12 - he12
end




function tprocedure!(conf::Config, grid::DeviceGrid{T1, T2}, mpts::DeviceParticle{T1, T2}) where {T1, T2}
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
        tresetgridstatus!(dev_grid)
        tresetmpstatus!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, conf.basis)
        tp2g!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts)
        tsolvegrid!(dev)(ndrange=dev_grid.ni, dev_grid, Δt)
        tdoublemapping1!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, Δt)
        tdoublemapping2!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts)
        tdoublemapping3!(dev)(ndrange=dev_grid.ni, dev_grid, Δt)
        tg2p!(dev)(ndrange=dev_mpts.np, dev_grid, dev_mpts, conf.material, t_eld, t_cur, Δt)
        t_cur += Δt
        h5.iters[] += 1
        update_pb!(printer, t_cur, t_tol)
    end
    finish_pb!(conf, printer); KAsync(dev)
    device2host!(mpts, dev_mpts)
    hdf5!(h5, fid, grid)
    close(fid)
end