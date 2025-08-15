export liE!
export dpP!
export material!

@inline Base.@propagate_inbounds function material!(
    mpts::DeviceParticle{T1, T2}, t_eld::T2, t_cur::T2, Δt::T2, 
    df1::T2, df2::T2, df3::T2, df4::T2, df5::T2, df6::T2, df7::T2, df8::T2, df9::T2,
    ix::Int, ::LinearElastic
) where {T1, T2}
    liE!(mpts, df1, df2, df3, df4, df5, df6, df7, df8, df9, ix)
end

@inline Base.@propagate_inbounds function material!(
    mpts::DeviceParticle{T1, T2}, t_eld::T2, t_cur::T2, Δt::T2, 
    df1::T2, df2::T2, df3::T2, df4::T2, df5::T2, df6::T2, df7::T2, df8::T2, df9::T2,
    ix::Int, ::DruckerPrager
) where {T1, T2}
    liE!(mpts, df1, df2, df3, df4, df5, df6, df7, df8, df9, ix)
    if t_cur > t_eld || t_eld == 0 
        dpP!(mpts, ix)
    end
end

@inline Base.@propagate_inbounds function liE!(mpts::DeviceParticle{T1, T2}, 
    df1::T2, df2::T2, df3::T2, df4::T2, df5::T2, df6::T2, df7::T2, df8::T2, df9::T2, ix::Int
) where {T1, T2}
    nid = mpts.nid[ix]
    Ks  = mpts.Ks[nid]
    Gs  = mpts.Gs[nid]
    # spin tensor
    # here, ωij = (vorticity tensor) × Δt, i.e.
    # ωij = (∂Ni × Vj - ∂Nj × Vi) × 0.5 × Δt
    ωxy = T2(0.5) * (df4 - df2)
    ωyz = T2(0.5) * (df8 - df6)
    ωxz = T2(0.5) * (df7 - df3)
    # objective stress
    # σij,new = σij,old + σij,R
    # σxx,R =  2 × (σxy × ωxy + σxz × ωxz)
    # σyy,R = -2 × (σxy × ωxy - σyz × ωyz)
    # σzz,R = -2 × (σxz × ωxz + σyz × ωyz)
    # σxy,R =  ωxy × (σyy - σxx) + σyz × ωxz + σxz × ωyz
    # σyz,R =  ωyz × (σzz - σyy) - σxy × ωxz - σxz × ωxy
    # σxz,R =  ωxz × (σzz - σxx) + σyz × ωxy - σxy × ωyz
    # where σij = σji, ωij = -ωji
    σij1 = mpts.σij[ix, 1]
    σij2 = mpts.σij[ix, 2]
    σij3 = mpts.σij[ix, 3]
    σij4 = mpts.σij[ix, 4]
    σij5 = mpts.σij[ix, 5]
    σij6 = mpts.σij[ix, 6]
    mpts.σij[ix, 1] +=  T2(2.0) * (σij4 * ωxy + σij6 * ωxz)
    mpts.σij[ix, 2] += -T2(2.0) * (σij4 * ωxy - σij5 * ωyz)
    mpts.σij[ix, 3] += -T2(2.0) * (σij6 * ωxz + σij5 * ωyz)
    mpts.σij[ix, 4] += ωxy * (σij2 - σij1) + ωxz * σij5 + ωyz * σij6
    mpts.σij[ix, 5] += ωyz * (σij3 - σij2) - ωxz * σij4 - ωxy * σij6
    mpts.σij[ix, 6] += ωxz * (σij3 - σij1) + ωxy * σij5 - ωyz * σij4
    # linear elastic
    Dt = Ks + T2(1.333333) * Gs
    Dd = Ks - T2(0.666667) * Gs
    mpts.σij[ix, 1] += Dt * df1 + Dd * df5 + Dd * df9
    mpts.σij[ix, 2] += Dd * df1 + Dt * df5 + Dd * df9
    mpts.σij[ix, 3] += Dd * df1 + Dd * df5 + Dt * df9
    mpts.σij[ix, 4] += Gs * (df2 + df4)
    mpts.σij[ix, 5] += Gs * (df6 + df8)
    mpts.σij[ix, 6] += Gs * (df3 + df7)
end

@inline Base.@propagate_inbounds function dpP!(mpts::DeviceParticle{T1, T2}, ix::Int) where {T1, T2}
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
        Δλs  = fs / (Gs + Ks * qϕ * qψ)
        tmp1 = σm - Ks * qψ * Δλs
        tmp2 = (kϕ - qϕ * tmp1) / τ
        mpts.σij[ix, 1] = s1 * tmp2 + tmp1
        mpts.σij[ix, 2] = s2 * tmp2 + tmp1
        mpts.σij[ix, 3] = s3 * tmp2 + tmp1
        mpts.σij[ix, 4] = s4 * tmp2
        mpts.σij[ix, 5] = s5 * tmp2
        mpts.σij[ix, 6] = s6 * tmp2
        mpts.ϵq[ix] += Δλs * sqrt(T2(0.333333) + T2(0.222222) * qψ * qψ)
        mpts.ϵk[ix] += Δλs * qψ
    end
    ## tensile failure correction
    if (σm ≥ σt) && (BF ≤ T2(0.0))
        Δλt = ft / Ks
        mpts.σij[ix, 1] = s1 + σt
        mpts.σij[ix, 2] = s2 + σt
        mpts.σij[ix, 3] = s3 + σt
        mpts.ϵq[ix] += Δλt * T2(0.333333) * T2(1.414214)
        mpts.ϵk[ix] += Δλt
    end
end