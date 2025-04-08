#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : timestep.jl                                                                |
|  Description: Adaptive time step algorithms in MaterialPointSolver.jl                    |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  License    : MIT License                                                                |
+==========================================================================================#

export Δtu, Δtl

function Δtu(Gs, Ks, Kw, ρs, ρw, n, dx)
    Ec = Ks + inv((n/Kw) + ((1-n)/Ks)) + 4/3*Gs
    ωu = 2/dx * sqrt((Ec + Kw/n)/((1-n)*ρs + n*ρw))
    return 2/ωu
end

function Δtl(Kw, k, ρw, n, μ, dx)
    ωl = 2/dx * sqrt(Kw/ρw)
    ξl = ((n*μ) / (2*k*ωl)) * inv(ρw)
    return 2/ωl * (sqrt(ξl^2 + 1) - ξl)
end