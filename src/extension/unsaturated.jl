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

export SWCC, HCC

@inline Base.@propagate_inbounds function SWCC(
    S_min::T2, 
    S_max::T2,
    P    ::T2, 
    P_ref::T2, 
    λ    ::T2
)::T2 where T2
    # ref [1]
    return S_min + (S_max - S_min) * (T2(1.0) + (P / P_ref)^inv(T2(1.0) - λ))^(-λ)
end

@inline Base.@propagate_inbounds function HCC(
    k_sat::T2, 
    α    ::T2, 
    n    ::T2, 
    P    ::T2
)::T2 where T2
    numtr = (T2(1.0) - (α*P)^(n - T2(1.0)) * (T2(1.0) + (α*P)^n)^(inv(n) - T2(1.0)))^2
    dentr = (T2(1.0) + (α*P)^n)^(T2(0.5) - T2(0.5) * n)
    
    # ref [2]
    return k_sat * (numtr / dentr)
end
