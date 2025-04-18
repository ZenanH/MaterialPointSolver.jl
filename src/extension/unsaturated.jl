#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : unsaturated.jl                                                             |
|  Description: SWCC and HCC curve functions                                               |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 07/10/2024                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  References : [1] Ceccato, F., Yerro, A., Girardi, V., Simonini, P., 2021. Two-phase     |
|                   dynamic MPM formulation for unsaturated soil. Computers and Geotechnics|
|                   129, 103876. https://doi.org/10.1016/j.compgeo.2020.103876             |
|               [2] Wang, D., Wang, B., Yuan, W., Liu, L., 2023. Investigation of rainfall |
|                   intensity on the slope failure process using GPU-accelerated coupled   |
|                   MPM. Comput. Geotech. 163, 105718.                                     |
|                   https://doi.org/10.1016/j.compgeo.2023.105718                          |
+==========================================================================================#

export SWCC, HCC, ∂SWCC

@inline Base.@propagate_inbounds function SWCC(
    S_min::T2, # min saturation
    S_max::T2, # max saturation
    P    ::T2, # pore pressure
    P_ref::T2, # fig parameter 1
    λ    ::T2  # fit parameter 2
)::T2 where T2
    p = P ≤ T2(0.0) ? T2(0.0) : P    
    # ref [1]
    return S_min + (S_max - S_min) * (T2(1.0) + (p / P_ref)^inv(T2(1.0) - λ))^(-λ)
end

@inline Base.@propagate_inbounds function HCC(
    k_sat::T2, 
    α    ::T2, 
    β    ::T2, 
    P    ::T2
)::T2 where T2
    p = P ≤ T2(0.0) ? T2(0.0) : P
    numtr = (T2(1.0) - (α*p)^(β - T2(1.0)) * (T2(1.0) + (α*p)^β)^(inv(β) - T2(1.0)))^2
    dentr = (T2(1.0) + (α*p)^β)^(T2(0.5) - T2(0.5) * β)
    
    # ref [2]
    return k_sat * (numtr / dentr)
end

@inline Base.@propagate_inbounds function ∂SWCC(
    S_min::T2, # min saturation
    S_max::T2, # max saturation
    P    ::T2, # pore pressure
    P_ref::T2, # fitting parameter 1
    λ    ::T2  # fitting parameter 2
)::T2 where T2
    p = P ≤ T2(0.0) ? T2(0.0) : P
    tmp1 = -λ * (S_max - S_min) * inv(T2(1.0) - λ)
    tmp2 = p^(λ*inv(T2(1.0) - λ)) / P_ref^inv(T2(1.0) - λ)
    tmp3 = (T2(1.0) + (p/P_ref)^inv(T2(1.0) - λ))^(-λ - T2(1.0))
    return tmp1 * tmp2 * tmp3
end