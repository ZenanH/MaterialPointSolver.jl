export bspline2basis
export linearbasis
export uGIMPbasis

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

@inline Base.@propagate_inbounds function linearbasis(x1::T2, x2::T2, invh::T2) where T2
    N1 = T2(1.0) - x1 * invh
    d1 = -invh
    N2 = T2(1.0) - x2 * invh
    d2 = invh
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