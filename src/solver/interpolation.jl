export bspline2basis

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