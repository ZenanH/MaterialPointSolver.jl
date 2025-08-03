#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: Calculator funcs                                                           |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, ISTE, Université de Lausanne                                   |
|  Maintainer : Zenan Huo                                                                  |
+==========================================================================================#

export update_F!, detF

@inline Base.@propagate_inbounds function update_F!(
    mp, df1, df2, df3, df4, T2_1, ix, dim::Dim2D
)
    F1 = mp.F[1, ix]; F2 = mp.F[2, ix]; F3 = mp.F[3, ix]; F4 = mp.F[4, ix]      
    mp.F[1, ix] = (df1 + T2_1) * F1 + df2 * F3
    mp.F[2, ix] = (df1 + T2_1) * F2 + df2 * F4
    mp.F[3, ix] = (df4 + T2_1) * F3 + df3 * F1
    mp.F[4, ix] = (df4 + T2_1) * F4 + df3 * F2
end

@inline Base.@propagate_inbounds function update_F!(
    mp, df1, df2, df3, df4, df5, df6, df7, df8, df9, T2_1, ix, dim::Dim3D
)
        F1 = mp.F[1, ix]; F2 = mp.F[2, ix]; F3 = mp.F[3, ix]
        F4 = mp.F[4, ix]; F5 = mp.F[5, ix]; F6 = mp.F[6, ix]
        F7 = mp.F[7, ix]; F8 = mp.F[8, ix]; F9 = mp.F[9, ix]        
        mp.F[1, ix] = (df1 + T2_1) * F1 + df2 * F4 + df3 * F7
        mp.F[2, ix] = (df1 + T2_1) * F2 + df2 * F5 + df3 * F8
        mp.F[3, ix] = (df1 + T2_1) * F3 + df2 * F6 + df3 * F9
        mp.F[4, ix] = (df5 + T2_1) * F4 + df4 * F1 + df6 * F7
        mp.F[5, ix] = (df5 + T2_1) * F5 + df4 * F2 + df6 * F8
        mp.F[6, ix] = (df5 + T2_1) * F6 + df4 * F3 + df6 * F9
        mp.F[7, ix] = (df9 + T2_1) * F7 + df8 * F4 + df7 * F1
        mp.F[8, ix] = (df9 + T2_1) * F8 + df8 * F5 + df7 * F2
        mp.F[9, ix] = (df9 + T2_1) * F9 + df8 * F6 + df7 * F3
end


@inline Base.@propagate_inbounds function detF(mp, ix, dim::Dim2D)
    return mp.F[1, ix] * mp.F[4, ix] - mp.F[2, ix] * mp.F[3, ix]
end

@inline Base.@propagate_inbounds function detF(mp, ix, dim::Dim3D)
    return mp.F[1, ix] * mp.F[5, ix] * mp.F[9, ix] + mp.F[2, ix] * mp.F[6, ix] * mp.F[7, ix] +
           mp.F[3, ix] * mp.F[4, ix] * mp.F[8, ix] - mp.F[7, ix] * mp.F[5, ix] * mp.F[3, ix] -
           mp.F[8, ix] * mp.F[6, ix] * mp.F[1, ix] - mp.F[9, ix] * mp.F[4, ix] * mp.F[2, ix]
end