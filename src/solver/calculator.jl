export update_F!
export detF

@inline Base.@propagate_inbounds function update_F!(mpts::DeviceParticle{T1, T2}, 
    df1::T2, df2::T2, df3::T2, df4::T2, df5::T2, df6::T2, df7::T2, df8::T2, df9::T2, ix::Int
) where {T1, T2}
    F1 = mpts.F[ix, 1]; F2 = mpts.F[ix, 2]; F3 = mpts.F[ix, 3]
    F4 = mpts.F[ix, 4]; F5 = mpts.F[ix, 5]; F6 = mpts.F[ix, 6]
    F7 = mpts.F[ix, 7]; F8 = mpts.F[ix, 8]; F9 = mpts.F[ix, 9]        
    mpts.F[ix, 1] = (df1 + T2(1.0)) * F1 + df2 * F4 + df3 * F7
    mpts.F[ix, 2] = (df1 + T2(1.0)) * F2 + df2 * F5 + df3 * F8
    mpts.F[ix, 3] = (df1 + T2(1.0)) * F3 + df2 * F6 + df3 * F9
    mpts.F[ix, 4] = (df5 + T2(1.0)) * F4 + df4 * F1 + df6 * F7
    mpts.F[ix, 5] = (df5 + T2(1.0)) * F5 + df4 * F2 + df6 * F8
    mpts.F[ix, 6] = (df5 + T2(1.0)) * F6 + df4 * F3 + df6 * F9
    mpts.F[ix, 7] = (df9 + T2(1.0)) * F7 + df8 * F4 + df7 * F1
    mpts.F[ix, 8] = (df9 + T2(1.0)) * F8 + df8 * F5 + df7 * F2
    mpts.F[ix, 9] = (df9 + T2(1.0)) * F9 + df8 * F6 + df7 * F3
end

@inline Base.@propagate_inbounds function detF(mpts::DeviceParticle{T1, T2}, ix::Int) where {T1, T2}
    return mpts.F[ix, 1] * mpts.F[ix, 5] * mpts.F[ix, 9] + mpts.F[ix, 2] * mpts.F[ix, 6] * mpts.F[ix, 7] +
           mpts.F[ix, 3] * mpts.F[ix, 4] * mpts.F[ix, 8] - mpts.F[ix, 7] * mpts.F[ix, 5] * mpts.F[ix, 3] -
           mpts.F[ix, 8] * mpts.F[ix, 6] * mpts.F[ix, 1] - mpts.F[ix, 9] * mpts.F[ix, 4] * mpts.F[ix, 2]
end