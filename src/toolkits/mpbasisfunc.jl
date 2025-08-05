#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : mpbasisfunc.jl                                                             |
|  Description: Basis functions for MPM                                                    |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Functions  : 1. linearbasis                                                             |
|               2. uGIMPbasis                                                              |
|               3. bspline2basis                                                           |
|               4. bspline3basis                                                           |
|               4. get_type                                                                |
+==========================================================================================#

export linearbasis, uGIMPbasis, bspline2basis, bspline3basis
export get_type

"""
    linearbasis(x1::T2, x2, h::T2)

Description:
---
Standard `linear basis function` for MPM, where `x1`, `x2` are the distance between particle
and node, and `h` is the grid spacing.
"""
@inline Base.@propagate_inbounds function linearbasis(x1::T2, x2::T2, h::T2) where T2
    h_denom = inv(h)
    N1 = T2(1.0) - x1 * h_denom
    d1 = -h_denom
    N2 = T2(1.0) - x2 * h_denom
    d2 = h_denom
    return N1, N2, d1, d2
end

@inline Base.@propagate_inbounds function uGIMPbasis(
    x1::T2, 
    x2::T2, 
    x3::T2, 
    h ::T2, 
    lp::T2
) where T2
    N1 = d1 = N2 = d2 = N3 = d3 = T2(0.0)
    absx1, absx2, absx3 = abs(x1), abs(x2), abs(x3)

    if absx1 < T2(0.5)*lp
        N1 = T2(1.0) - ((T2(4.0) * x1 * x1 + lp * lp) / (T2(4.0) * h * lp))
        d1 = -((T2(8.0) * x1) / (T2(4.0) * h * lp))
    elseif (T2(0.5) * lp) ≤ absx1 < (h - T2(0.5) * lp)
        N1 = T2(1.0) - (absx1 / h)
        d1 = sign(x1) * (T2(-1.0) / h)
    elseif (h - T2(0.5) * lp) ≤ absx1 < (h + T2(0.5) * lp)
        N1 = ((h + T2(0.5) * lp - absx1) * (h + T2(0.5) * lp - absx1)) / (T2(2.0) * h * lp)
        d1 = -sign(x1) * ((h + T2(0.5) * lp - absx1) / (h * lp))
    end

    if absx2 < T2(0.5)*lp
        N2 = T2(1.0) - ((T2(4.0) * x2 * x2 + lp * lp) / (T2(4.0) * h * lp))
        d2 = -((T2(8.0) * x2) / (T2(4.0) * h * lp))
    elseif (T2(0.5) * lp) ≤ absx2 < (h - T2(0.5) * lp)
        N2 = T2(1.0) - (absx2 / h)
        d2 = sign(x2) * (T2(-1.0) / h)
    elseif (h - T2(0.5) * lp) ≤ absx2 < (h + T2(0.5) * lp)
        N2 = ((h + T2(0.5) * lp - absx2) * (h + T2(0.5) * lp - absx2)) / (T2(2.0) * h * lp)
        d2 = -sign(x2) * ((h + T2(0.5) * lp - absx2) / (h * lp))
    end

    if absx3 < T2(0.5)*lp
        N3 = T2(1.0) - ((T2(4.0) * x3 * x3 + lp * lp) / (T2(4.0) * h * lp))
        d3 = -((T2(8.0) * x3) / (T2(4.0) * h * lp))
    elseif (T2(0.5) * lp) ≤ absx3 < (h - T2(0.5) * lp)
        N3 = T2(1.0) - (absx3 / h)
        d3 = sign(x3) * (T2(-1.0) / h)
    elseif (h - T2(0.5) * lp) ≤ absx3 < (h + T2(0.5) * lp)
        N3 = ((h + T2(0.5) * lp - absx3) * (h + T2(0.5) * lp - absx3)) / (T2(2.0) * h * lp)
        d3 = -sign(x3) * ((h + T2(0.5) * lp - absx3) / (h * lp))
    end
    
    return N1, N2, N3, d1, d2, d3
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
    return T2(Ni), T2(dN * inv(h))
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