#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : bingham.jl                                                                 |
|  Description: Bingham plastic material's constitutive model                              |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Functions  : 1. bhP! [2D]                                                               |
|               2. bhP! [3D]                                                               |
+==========================================================================================#

export bhP!

@kernel inbounds = true function bhP!(
    mp  ::DeviceParticle2D{T1, T2},
    attr::  DeviceProperty{T1, T2},
    ΔT_1::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        nid  = attr.nid[ix]
        Ks   = attr.Ks[nid]
        Gs   = attr.Gs[nid]
        μd   = attr.μd[nid]
        τy   = attr.τy[nid]
        ηmax = attr.ηmax[nid]
        # get strain rate for xx, yy, xy
        dϵxx = ΔT_1 * mp.ΔFs[ix, 1] 
        dϵyy = ΔT_1 * mp.ΔFs[ix, 4]
        dϵxy = ΔT_1 * (mp.ΔFs[ix, 2] + mp.ΔFs[ix, 3]) * T2(0.5)
        # compute elastic increments
        p = mp.σm[ix] + Ks * (mp.ΔFs[ix, 1] + mp.ΔFs[ix, 4])
        # elastic trial stress
        γ̇ = sqrt(T2(2.0) * (dϵxx * dϵxx + dϵyy * dϵyy + T2(2.0) * dϵxy * dϵxy))
        τxx = Gs * (dϵxx - dϵyy)
        τyy = Gs * (dϵyy - dϵxx)
        τxy = Gs *  dϵxy * T2(2.0)
        τc  = sqrt(T2(0.5) * (τxx * τxx + τyy * τyy + T2(2.0) * τxy * τxy))
        # bingham law
        ηp = T2(0.0)
        if τc > τy
            ηp = μd + τy / γ̇
            ηp = ηp > ηmax ? ηmax : ηp
        end
        sxx = T2(2.0) * ηp * dϵxx
        syy = T2(2.0) * ηp * dϵyy
        sxy = T2(2.0) * ηp * dϵxy
        # update stress
        mp.σij[ix, 1] = sxx + p
        mp.σij[ix, 2] = syy + p
        mp.σij[ix, 4] = sxy
        mp.σm[ix] = (mp.σij[ix, 1] + mp.σij[ix, 2]) * T2(0.5)
        mp.sij[ix, 1] = mp.σij[ix, 1] - mp.σm[ix]
        mp.sij[ix, 2] = mp.σij[ix, 2] - mp.σm[ix]
        mp.sij[ix, 4] = mp.σij[ix, 4]
    end
end

@kernel inbounds = true function bhP!(
    mp  ::DeviceParticle3D{T1, T2},
    attr::  DeviceProperty{T1, T2},
    ΔT_1::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        nid  = attr.nid[ix]
        Ks   = attr.Ks[nid]
        Gs   = attr.Gs[nid]
        μd   = attr.μd[nid]
        τy   = attr.τy[nid]
        ηmax = attr.ηmax[nid]
        # get strain rate for xx, yy, zz, xy, yz, zx
        dϵxx = ΔT_1 * mp.ΔFs[ix, 1] 
        dϵyy = ΔT_1 * mp.ΔFs[ix, 5]
        dϵzz = ΔT_1 * mp.ΔFs[ix, 9]
        dϵxy = ΔT_1 * (mp.ΔFs[ix, 2] + mp.ΔFs[ix, 4]) * T2(0.5)
        dϵyz = ΔT_1 * (mp.ΔFs[ix, 6] + mp.ΔFs[ix, 8]) * T2(0.5)
        dϵzx = ΔT_1 * (mp.ΔFs[ix, 3] + mp.ΔFs[ix, 7]) * T2(0.5)
        # compute elastic increments
        p = mp.σm[ix] + Ks * (mp.ΔFs[ix, 1] + mp.ΔFs[ix, 5] + mp.ΔFs[ix, 9])
        # elastic trial stress
        γ̇ = sqrt(T2(2.0) * (dϵxx * dϵxx + dϵyy * dϵyy + dϵzz * dϵzz + 
            T2(2.0) * dϵxy * dϵxy + T2(2.0) * dϵyz * dϵyz + T2(2.0) * dϵzx * dϵzx))
        trdϵ = T2(0.333333) * (dϵxx + dϵyy + dϵzz)
        τxx = T2(2.0) * Gs * (dϵxx - trdϵ)
        τyy = T2(2.0) * Gs * (dϵyy - trdϵ)
        τzz = T2(2.0) * Gs * (dϵzz - trdϵ)
        τxy = T2(2.0) * Gs * dϵxy
        τyz = T2(2.0) * Gs * dϵyz
        τzx = T2(2.0) * Gs * dϵzx
        τc  = sqrt(T2(0.5) * ((τxx - τyy) * (τxx - τyy) + 
                              (τyy - τzz) * (τyy - τzz) + 
                              (τzz - τxx) * (τzz - τxx) + T2(6.0) * 
            (τxy * τxy + τyz * τyz + τzx * τzx)))
        # bingham law
        ηp = T2(0.0)
        if τc > τy
            ηp = μd + τy / γ̇
            ηp = ηp > ηmax ? ηmax : ηp
        end
        sxx = T2(2.0) * ηp * dϵxx
        syy = T2(2.0) * ηp * dϵyy
        szz = T2(2.0) * ηp * dϵzz
        sxy = T2(2.0) * ηp * dϵxy
        syz = T2(2.0) * ηp * dϵyz
        szx = T2(2.0) * ηp * dϵzx
        # update stress
        mp.σij[ix, 1] = sxx + p
        mp.σij[ix, 2] = syy + p
        mp.σij[ix, 3] = szz + p
        mp.σij[ix, 4] = sxy
        mp.σij[ix, 5] = syz
        mp.σij[ix, 6] = szx
        mp.σm[ix] = (mp.σij[ix, 1] + mp.σij[ix, 2] + mp.σij[ix, 3]) * T2(0.333333)
        mp.sij[ix, 1] = mp.σij[ix, 1] - mp.σm[ix]
        mp.sij[ix, 2] = mp.σij[ix, 2] - mp.σm[ix]
        mp.sij[ix, 3] = mp.σij[ix, 3] - mp.σm[ix]
        mp.sij[ix, 4] = mp.σij[ix, 4]
        mp.sij[ix, 5] = mp.σij[ix, 5]
        mp.sij[ix, 6] = mp.σij[ix, 6]
    end
end

# # Herschel-Bulkley model for 2D
# @kernel inbounds = true function testbhP2!(
#     mp  ::DeviceParticle2D{T1, T2},
#     attr::  DeviceProperty{T1, T2},
#     ΔT_1::T2
# ) where {T1, T2}
#     ix = @index(Global)
#     if ix ≤ mp.np
#         nid  = attr.nid[ix]
#         Ks   = attr.Ks[nid]
#         Gs   = attr.Gs[nid]
#         μd   = attr.μd[nid]
#         τy   = attr.τy[nid]
#         ηmax = attr.ηmax[nid]
#         # get strain rate for xx, yy, xy
#         dϵxx = ΔT_1 * mp.ΔFs[ix, 1] 
#         dϵyy = ΔT_1 * mp.ΔFs[ix, 4]
#         dϵxy = ΔT_1 * (mp.ΔFs[ix, 2] + mp.ΔFs[ix, 3]) * T2(0.5)
#         # compute elastic increments
#         p = mp.σm[ix] + Ks * (mp.ΔFs[ix, 1] + mp.ΔFs[ix, 4])
#         # elastic trial stress
#         γ̇ = sqrt(T2(2.0) * (dϵxx * dϵxx + dϵyy * dϵyy + T2(2.0) * dϵxy * dϵxy))
#         τxx = Gs * (dϵxx - dϵyy)
#         τyy = Gs * (dϵyy - dϵxx)
#         τxy = Gs *  dϵxy * T2(2.0)
#         τc  = sqrt(T2(0.5) * (τxx * τxx + τyy * τyy + T2(2.0) * τxy * τxy))
#         # bingham law
#         ηp = T2(0.0)
#         if τc > τy
#             K = 0.99
#             n = 1
#             ηp = (τy / γ̇) + K * γ̇^(n-1)
#             ηp = ηp > ηmax ? ηmax : ηp
#         end
#         sxx = T2(2.0) * ηp * dϵxx
#         syy = T2(2.0) * ηp * dϵyy
#         sxy = T2(2.0) * ηp * dϵxy
#         # update stress
#         mp.σij[ix, 1] = sxx + p
#         mp.σij[ix, 2] = syy + p
#         mp.σij[ix, 4] = sxy
#         mp.σm[ix] = (mp.σij[ix, 1] + mp.σij[ix, 2]) * T2(0.5)
#     end
# end