#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : timestep.jl                                                                |
|  Description: Adaptive time step algorithms in MaterialPointSolver.jl                    |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  License    : MIT License                                                                |
|  References : [1]. Colom, A.Y., n.d. MPM modelling of landslides in brittle and          |
|               unsaturated soils.                                                         |
|               [2]. Yerro, A., Girardi, V., Martinelli, M., Ceccato, F., 2022. Modelling  |
|               unsaturated soils with the material point method. A discussion of the      |
|               state-of-the-art. Geomechanics for Energy and the Environment 32, 100343.  |
|               https://doi.org/10.1016/j.gete.2022.100343.                                |
+==========================================================================================#

export satΔt, unsatΔt

@inline function satΔt(
    Gs::T2, 
    Ks::T2, 
    Kw::T2, 
    ρs::T2, 
    ρw::T2, 
    n ::T2, 
    k ::T2, 
    h ::T2
) where T2
    # sat_Δt,u
    Ec = Ks + 1/((n/Kw) + ((1-n)/Ks)) + (4/3)*Gs
    ωu = 2/h * sqrt((Ec + Kw/n)/((1-n)*ρs + n*ρw))
    sat_Δtu = 2/ωu

    # sat_Δt,l
    ωl = 2/h * sqrt(Kw/ρw)
    ξl = ((n*ρw*T2(9.8)) / (2*k*ωl)) * (1/ρw)
    sat_Δtl = 2/ωl * (sqrt(ξl^2 + 1) - ξl)

    # ref[1], eq. 4.84
    return T2(min(sat_Δtu, sat_Δtl))
end

@inline function unsatΔt(
    Gs::T2, 
    Ks::T2, 
    Kw::T2, 
    ρs::T2, 
    ρw::T2, 
    n ::T2, 
    k ::T2, 
    h ::T2, 
    S ::T2, 
    λ ::T2
) where T2
    # Δtc Mieremet
    λc = S/(λ + S/Kw) # λ=∂S/∂s (S: degree of saturation, s: suction ↣ ∂S/∂p (pore pressure))
    ρm = (1-n)*ρs + n*S*ρw
    Ec = Ks + 4/3 * Gs

    # parameters
    a = (S^2 * n * ρm * ρw * T2(9.8)) / ((1-n) * ρs * ρw * k)
    b = (4*(n*ρm*λc + (1-2n)*ρw*λc + n*ρw*Ec)) / (n*(1-n)*ρs*ρw*(h/2)^2)
    d = (16*Ec*λc) / ((1-n)*ρs*ρw*(h/2)^4)

    # ref[2], eq. 35
    return T2((-2a + sqrt(4a^2 + 8(b + sqrt(b^2 - 4d)))) / (b + sqrt(b^2 - 4d)))
end