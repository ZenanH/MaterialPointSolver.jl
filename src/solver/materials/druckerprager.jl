#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: Drucker-Prager material's constitutive model.                              |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, ISTE, Université de Lausanne                                   |
|  Maintainer : Zenan Huo                                                                  |
+==========================================================================================#

export dpP!

@inline Base.@propagate_inbounds function dpP!(
    mpts, ix, dim::Dim2D, ϵ::Precision{T1, T2}
) where {T1, T2}
    nid = mpts.nid[ix]
    c   = mpts.c[nid]
    ϕ   = mpts.ϕ[nid]
    ψ   = mpts.ψ[nid]
    σt  = mpts.σt[nid]
    Gs  = mpts.Gs[nid]
    Ks  = mpts.Ks[nid]
    # update deviatoric stress tensor
    σm = (mpts.σij[1, ix] + mpts.σij[2, ix] + mpts.σij[3, ix]) * T2(0.333333)
    s1 = mpts.σij[1, ix] - σm
    s2 = mpts.σij[2, ix] - σm
    s3 = mpts.σij[3, ix] - σm
    s4 = mpts.σij[4, ix]

    # drucker-prager
    τ = sqrt(T2(0.5) * (s1 * s1 + s2 * s2 + s3 * s3) + s4 * s4)
    kϕ = (T2(6.0) * c * cos(ϕ)) / (T2(1.732051) * (T2(3.0) + sin(ϕ)))
    qϕ = (T2(6.0)     * sin(ϕ)) / (T2(1.732051) * (T2(3.0) + sin(ϕ)))
    qψ = (T2(6.0)     * sin(ψ)) / (T2(1.732051) * (T2(3.0) + sin(ψ)))
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
        mpts.σij[1, ix] = s1 * tmp2 + tmp1
        mpts.σij[2, ix] = s2 * tmp2 + tmp1
        mpts.σij[3, ix] = s3 * tmp2 + tmp1
        mpts.σij[4, ix] = s4 * tmp2
        mpts.ϵq[ix] += Δλs * sqrt(T2(0.333333) + T2(0.222222) * qψ *qψ)
        mpts.ϵk[ix] += Δλs * qψ
    end
    ## tensile failure correction
    if (σm ≥ σt) && (BF ≤ T2(0.0))
        Δλt = ft / Ks
        mpts.σij[1, ix] = s1 + σt
        mpts.σij[2, ix] = s2 + σt
        mpts.σij[3, ix] = s3 + σt
        mpts.ϵq[ix] += Δλt * T2(0.333333) * T2(1.414214)
        mpts.ϵk[ix] += Δλt
    end
end

@inline Base.@propagate_inbounds function dpP!(
    mpts, ix, dim::Dim3D, ϵ::Precision{T1, T2}
) where {T1, T2}
    nid = mpts.nid[ix]
    c   = mpts.c[nid]
    ϕ   = mpts.ϕ[nid]
    ψ   = mpts.ψ[nid]
    σt  = mpts.σt[nid]
    Gs  = mpts.Gs[nid]
    Ks  = mpts.Ks[nid]
    # update deviatoric stress tensor
    σm = (mpts.σij[1, ix] + mpts.σij[2, ix] + mpts.σij[3, ix]) * T2(0.333333)
    s1 = mpts.σij[1, ix] - σm
    s2 = mpts.σij[2, ix] - σm
    s3 = mpts.σij[3, ix] - σm
    s4 = mpts.σij[4, ix]
    s5 = mpts.σij[5, ix]
    s6 = mpts.σij[6, ix]

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
        mpts.σij[1, ix] = s1 * tmp2 + tmp1
        mpts.σij[2, ix] = s2 * tmp2 + tmp1
        mpts.σij[3, ix] = s3 * tmp2 + tmp1
        mpts.σij[4, ix] = s4 * tmp2
        mpts.σij[5, ix] = s5 * tmp2
        mpts.σij[6, ix] = s6 * tmp2
        mpts.ϵq[ix] += Δλs * sqrt(T2(0.333333) + T2(0.222222) * qψ * qψ)
        mpts.ϵk[ix] += Δλs * qψ
    end
    ## tensile failure correction
    if (σm ≥ σt) && (BF ≤ T2(0.0))
        Δλt = ft / Ks
        mpts.σij[1, ix] = s1 + σt
        mpts.σij[2, ix] = s2 + σt
        mpts.σij[3, ix] = s3 + σt
        mpts.ϵq[ix] += Δλt * T2(0.333333) * T2(1.414214)
        mpts.ϵk[ix] += Δλt
    end
end