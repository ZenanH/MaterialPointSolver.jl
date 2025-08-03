#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: Basis (shape) funcs                                                        |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, ISTE, Université de Lausanne                                   |
|  Maintainer : Zenan Huo                                                                  |
+==========================================================================================#

export linearbasis

@inline Base.@propagate_inbounds function linearbasis(x1::T2, x2::T2, invh::T2) where T2
    N1 = T2(1.0) - x1 * invh
    d1 = -invh
    N2 = T2(1.0) - x2 * invh
    d2 = invh
    return N1, N2, d1, d2
end