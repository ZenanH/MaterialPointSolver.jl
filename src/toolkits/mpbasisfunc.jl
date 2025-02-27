#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : mpbasisfunc.jl                                                             |
|  Description: Basis functions for MPM                                                    |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Functions  : 1. linearBasis                                                             |
|               2. uGIMPBasis                                                              |
|               3. bspline2basis                                                           |
|               4. bspline3basis                                                           |
|               4. get_type                                                                |
+==========================================================================================#

export linearbasis, uGIMPbasis, bspline2basis, bspline3basis
export uGIMPbasisx, uGIMPbasisy, uGIMPbasisz
export get_type

"""
    linearbasis(Δx::T2, h::T2)

Description:
---
Standard `linear basis function` for MPM, where `Δx` is the distance between particle and 
    node, and `h` is the grid spacing.
"""
@inline Base.@propagate_inbounds function linearbasis(Δx::T2, h::T2) where T2
    h_denom = inv(h)
    Ni = T2(1.0) - abs(Δx) * h_denom
    dN = -sign(Δx) * h_denom
    return T2(Ni), T2(dN)
end

"""
    uGIMPbasis(Δx::T2, h::T2, lp::T2)

Description:
---
`Uniform Generalized Interpolation` MPM basis function for MPM, where `Δx` is the distance 
between particle and node, `h` is the grid spacing, `lp` is the particle spacing.
"""
# This version will take a longer time since too many registers are used.
@inline Base.@propagate_inbounds function uGIMPbasis(Δx::T2, h::T2, lp::T2) where T2
    Ni = dN = T2(0.0)
    if abs(Δx) < T2(0.5)*lp
        Ni = T2(1.0) - ((T2(4.0) * Δx * Δx + lp * lp) / (T2(4.0) * h * lp))
        dN = -((T2(8.0) * Δx) / (T2(4.0) * h * lp))
    elseif (T2(0.5) * lp) ≤ abs(Δx) < (h - T2(0.5) * lp)
        Ni = T2(1.0) - (abs(Δx) / h)
        dN = sign(Δx) * (T2(-1.0) / h)
    elseif (h - T2(0.5) * lp) ≤ abs(Δx) < (h + T2(0.5) * lp)
        Ni = ((h + T2(0.5) * lp - abs(Δx)) * (h + T2(0.5) * lp - abs(Δx))) / (T2(2.0) * h * lp)
        dN = -sign(Δx) * ((h + T2(0.5) * lp - abs(Δx)) / (h * lp))
    end
    return T2(Ni), T2(dN)
end

@inline Base.@propagate_inbounds function uGIMPbasisx(Δx::T2, smem) where T2
    # smem[1] = h
    # smem[2] = lp
    # smem[3] = T2(1.0) / (T2(4.0) * h * lp)
    # smem[4] = T2(1.0) / (h * lp)
    # smem[5] = T2(1.0) / h
    # smem[6] = T2(0.5) * lp
    absΔx = abs(Δx)
    Ni = dN = T2(0.0)
    if absΔx < smem[6]
        Ni = T2(1.0) - ((T2(4.0) * Δx * Δx + smem[2] * smem[2]) * smem[3])
        dN = -T2(8.0) * Δx * smem[3]
    elseif smem[6] ≤ absΔx < smem[1] - smem[6]
        Ni = T2(1.0) - (absΔx / smem[1])
        dN = sign(Δx) * -smem[5]
    elseif smem[1] - smem[6] ≤ absΔx < smem[1] + smem[6]
        Ni = (smem[1] + smem[6] - absΔx) * (smem[1] + smem[6] - absΔx) * T2(0.5) * smem[4]
        dN = -sign(Δx) * (smem[1] + smem[6] - absΔx) * smem[4]
    end
    return T2(Ni), T2(dN)
end

@inline Base.@propagate_inbounds function uGIMPbasisy(Δx::T2, smem) where T2
    # smem[7]  = h
    # smem[8]  = lp
    # smem[9]  = T2(1.0) / (T2(4.0) * h * lp)
    # smem[10] = T2(1.0) / (h * lp)
    # smem[11] = T2(1.0) / h
    # smem[12] = T2(0.5) * lp
    absΔx = abs(Δx)
    Ni = dN = T2(0.0)
    if absΔx < smem[12]
        Ni = T2(1.0) - ((T2(4.0) * Δx * Δx + smem[8] * smem[8]) * smem[9])
        dN = -T2(8.0) * Δx * smem[9]
    elseif smem[12] ≤ absΔx < smem[7] - smem[12]
        Ni = T2(1.0) - (absΔx / smem[7])
        dN = sign(Δx) * -smem[11]
    elseif smem[7] - smem[12] ≤ absΔx < smem[7] + smem[12]
        Ni = (smem[7] + smem[12] - absΔx) * (smem[7] + smem[12] - absΔx) * T2(0.5) * smem[10]
        dN = -sign(Δx) * (smem[7] + smem[12] - absΔx) * smem[10]
    end
    return T2(Ni), T2(dN)
end

@inline Base.@propagate_inbounds function uGIMPbasisz(Δx::T2, smem) where T2
    # smem[13] = h
    # smem[14] = lp
    # smem[15] = T2(1.0) / (T2(4.0) * h * lp)
    # smem[16] = T2(1.0) / (h * lp)
    # smem[17] = T2(1.0) / h
    # smem[18] = T2(0.5) * lp
    absΔx = abs(Δx)
    Ni = dN = T2(0.0)
    if absΔx < smem[18]
        Ni = T2(1.0) - ((T2(4.0) * Δx * Δx + smem[14] * smem[14]) * smem[15])
        dN = -T2(8.0) * Δx * smem[15]
    elseif smem[18] ≤ absΔx < smem[13] - smem[18]
        Ni = T2(1.0) - (absΔx / smem[13])
        dN = sign(Δx) * -smem[17]
    elseif smem[13] - smem[18] ≤ absΔx < smem[13] + smem[18]
        Ni = (smem[13] + smem[18] - absΔx) * (smem[13] + smem[18] - absΔx) * T2(0.5) * smem[16]
        dN = -sign(Δx) * (smem[13] + smem[18] - absΔx) * smem[16]
    end
    return T2(Ni), T2(dN)
end

@inline Base.@propagate_inbounds function bspline2basis(ξ1::T2, ξ2::T2, ξ3::T2) where T2
    N1 = T2(0.50) * (T2(1.5) - ξ1) * (T2(1.5) - ξ1)
    N2 = T2(0.75) - ξ2 * ξ2
    N3 = T2(0.50) * (T2(1.5) + ξ3) * (T2(1.5) + ξ3)
    return N1, N2, N3
end

@inline Base.@propagate_inbounds function bspline2basis(
    ξ1::T2,
    ξ2::T2, 
    ξ3::T2,
    Δ1::T2
) where T2
    N1 = T2(0.50) * (T2(1.5) - ξ1) * (T2(1.5) - ξ1)
    N2 = T2(0.75) - ξ2 * ξ2
    N3 = T2(0.50) * (T2(1.5) + ξ3) * (T2(1.5) + ξ3)
    d1 = ξ1 - T2(1.5)
    d2 = T2(-2.0) * ξ2
    d3 = ξ3 + T2(1.5)
    return N1, N2, N3, d1 * Δ1, d2 * Δ1, d3 * Δ1
end

@inline Base.@propagate_inbounds function bspline3basis(
    r     ::T2, 
    h     ::T2, 
    type_v::Int32
) where T2
    # 1/6 = 0.166667
    # 4/3 = 1.333333
    # 2/3 = 0.666667
    # 1/3 = 0.333333
    Ni, dN = T2(0.0), T2(0.0)
    if type_v == Int32(1)
        if Int32(-2) < r < Int32(-1)
            Ni = T2(0.166667) * r * r * r + r * r + r + r + T2(1.333333)
            dN = T2(0.5) * r * r + r + r + T2(2.0)
        elseif Int32(-1) ≤ r < Int32(0)
            Ni = T2(-0.166667) * r * r * r + r + T2(1.0)
            dN = T2(-0.5) * r * r + T2(1.0)
        elseif Int32(0) ≤ r < Int32(1)
            Ni = T2(0.166667) * r * r * r - r + T2(1.0)
            dN = T2(0.5) * r * r - T2(1.0)
        elseif Int32(1) ≤ r < Int32(2)
            Ni = T2(-0.166667) * r * r * r + r * r - r - r + T2(1.333333)
            dN = -T2(0.5) * r * r + r + r - T2(2.0)
        end
    elseif type_v == Int32(2)
        if Int32(-1) ≤ r < Int32(0)
            Ni = T2(-0.333333) * r * r * r - r * r + T2(0.666667)
            dN = - r * r - r - r
        elseif Int32(0) ≤ r < Int32(1)
            Ni = T2(0.5) * r * r * r - r * r + T2(0.666667)
            dN = T2(1.5) * r * r - r - r
        elseif Int32(1) ≤ r < Int32(2)
            Ni = T2(-0.166667) * r * r * r + r * r - r - r + T2(1.333333)
            dN = T2(-0.5) * r * r + r + r - T2(2.0)
        end
    elseif type_v == Int32(3)
        if Int32(-2) < r < Int32(-1)
            Ni = T2(0.166667) * r * r * r + r * r + r + r + T2(1.333333)
            dN = T2(0.5) * r * r + r + r + T2(2.0)
        elseif Int32(-1) ≤ r < Int32(0)
            Ni = T2(-0.5) * r * r * r - r * r + T2(0.666667)
            dN = T2(-1.5) * r * r - r - r
        elseif Int32(0) ≤ r < Int32(1)
            Ni = T2(0.5) * r * r * r - r * r + T2(0.666667)
            dN = T2(1.5) * r * r - r - r
        elseif Int32(1) ≤ r < Int32(2)
            Ni = T2(-0.166667) * r * r * r + r * r - r - r + T2(1.333333)
            dN = T2(-0.5) * r * r + r + r - T2(2.0)
        end 
    elseif type_v == Int32(4)
        if Int32(-2) < r < Int32(-1)
            Ni = T2(0.166667) * r * r * r + r * r + r + r + T2(1.333333)
            dN = T2(0.5) * r * r + r + r + T2(2.0)
        elseif Int32(-1) ≤ r < Int32(0)
            Ni = T2(-0.5) * r * r * r - r * r + T2(0.666667)
            dN = T2(-1.5) * r * r - r - r
        elseif Int32(0) ≤ r ≤ Int32(1)
            Ni = T2(0.333333) * r * r * r - r * r + T2(0.666667)
            dN = r * r - r - r
        end
    end
    return T2(Ni), T2(dN / h)
end

@inline Base.@propagate_inbounds function get_type(
    gξ::T2, 
    b1::T2, 
    b2::T2, 
    h ::T2
)::Int32 where T2
    if gξ == b1 || gξ == b2
        return Int32(1)
    elseif gξ == b1 + h
        return Int32(2)
    elseif gξ == b2 - h
        return Int32(4)
    else
        return Int32(3)
    end
end