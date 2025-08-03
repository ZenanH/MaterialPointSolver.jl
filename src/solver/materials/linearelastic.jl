#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: Linear elastic material's constitutive model.                              |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, ISTE, Université de Lausanne                                   |
|  Maintainer : Zenan Huo                                                                  |
+==========================================================================================#

export liE!

@inline Base.@propagate_inbounds function liE!(
    mpts, df1, df2, df3, df4, ix, dim::Dim2D, ϵ::Precision{T1, T2}
) where {T1, T2}
    nid = mpts.nid[ix]
    Ks  = mpts.Ks[nid]
    Gs  = mpts.Gs[nid]
    # spin tensor
    # here, ωij = (vorticity tensor) × Δt, i.e.
    # ωij = (∂Ni × Vj - ∂Nj × Vi) × 0.5 × Δt
    ωxy = T2(0.5) * (df2 - df3)
    # objective stress
    # σij,new = σij,old + σij,R
    # σxx,R =  2 × σxy × ωxy
    # σyy,R = -2 × σxy × ωxy
    # σxy,R =  ωxy × (σyy - σxx)
    # where σij = σji, ωij = -ωji
    σij1 = mpts.σij[1, ix]
    σij2 = mpts.σij[2, ix]
    σij4 = mpts.σij[4, ix]
    mpts.σij[1, ix] +=  ωxy *  σij4 * T2(2.0)
    mpts.σij[2, ix] += -ωxy *  σij4 * T2(2.0)
    mpts.σij[4, ix] +=  ωxy * (σij2 - σij1)
    # linear elastic
    Dt = Ks + T2(1.333333) * Gs
    Dd = Ks - T2(0.666667) * Gs
    mpts.σij[1, ix] += Dt * df1 + Dd * df4 + Dd * T2(0.0)
    mpts.σij[2, ix] += Dd * df1 + Dt * df4 + Dd * T2(0.0)
    mpts.σij[3, ix] += Dd * df1 + Dd * df4 + Dt * T2(0.0)
    mpts.σij[4, ix] += Gs * (df2 + df3)
end

@inline Base.@propagate_inbounds function liE!(
    mpts, df1, df2, df3, df4, df5, df6, df7, df8, df9, ix, dim::Dim3D, ϵ::Precision{T1, T2}
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
    σij1 = mpts.σij[1, ix]
    σij2 = mpts.σij[2, ix]
    σij3 = mpts.σij[3, ix]
    σij4 = mpts.σij[4, ix]
    σij5 = mpts.σij[5, ix]
    σij6 = mpts.σij[6, ix]
    mpts.σij[1, ix] +=  T2(2.0) * (σij4 * ωxy + σij6 * ωxz)
    mpts.σij[2, ix] += -T2(2.0) * (σij4 * ωxy - σij5 * ωyz)
    mpts.σij[3, ix] += -T2(2.0) * (σij6 * ωxz + σij5 * ωyz)
    mpts.σij[4, ix] += ωxy * (σij2 - σij1) + ωxz * σij5 + ωyz * σij6
    mpts.σij[5, ix] += ωyz * (σij3 - σij2) - ωxz * σij4 - ωxy * σij6
    mpts.σij[6, ix] += ωxz * (σij3 - σij1) + ωxy * σij5 - ωyz * σij4
    # linear elastic
    Dt = Ks + T2(1.333333) * Gs
    Dd = Ks - T2(0.666667) * Gs
    mpts.σij[1, ix] += Dt * df1 + Dd * df5 + Dd * df9
    mpts.σij[2, ix] += Dd * df1 + Dt * df5 + Dd * df9
    mpts.σij[3, ix] += Dd * df1 + Dd * df5 + Dt * df9
    mpts.σij[4, ix] += Gs * (df2 + df4)
    mpts.σij[5, ix] += Gs * (df6 + df8)
    mpts.σij[6, ix] += Gs * (df3 + df7)
end