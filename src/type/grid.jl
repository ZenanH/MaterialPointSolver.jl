export Grid, DeviceGrid
export init_grid

abstract type DeviceGrid{T1, T2} end

struct Grid{T1, T2,
    T3 <: AbstractArray{T1, 1},
    #T4 <: AbstractArray{T1, 2},
    T5 <: AbstractArray{T2, 1},
    T6 <: AbstractArray{T2, 2},
    T7 <: NamedTuple
} <: DeviceGrid{T1, T2}
    ni  ::T1
    nx  ::T1
    ny  ::T1
    nz  ::T1
    nc  ::T1
    ncx ::T1
    ncy ::T1
    ncz ::T1
    h   ::T2
    invh::T2
    x1  ::T2
    x2  ::T2
    y1  ::T2
    y2  ::T2
    z1  ::T2
    z2  ::T2
    ζs  ::T2
    ms  ::T5
    ps  ::T6
    vs  ::T6
    vsT ::T6
    fs  ::T6
    vsxi::T3
    vsxv::T5
    vsyi::T5
    vsyv::T5
    vszi::T5
    vszv::T5
    ext ::T7
end

@KAadapt Grid

function init_grid(xr::AbstractRange, yr::AbstractRange, zr::AbstractRange; 
    ϵ::Symbol=:double, vsxi=[1], vsxv=[0], vsyi=[1], vsyv=[0], vszi=[1], vszv=[0], 
    ζs=0.0, ext=[0]
)
    T1  = ϵ == :single ? Int32   : Int64
    T2  = ϵ == :single ? Float32 : Float64
    x, y, z = convert.(Float64, xr), convert.(Float64, yr), convert.(Float64, zr)
    nx = T1(length(x)); nx ≤ 1 && error("Grid must have at least 2 points in x-direction")
    ny = T1(length(y)); ny ≤ 1 && error("Grid must have at least 2 points in y-direction")
    nz = T1(length(z)); nz ≤ 1 && error("Grid must have at least 2 points in z-direction")
    ni = T1(nx * ny * nz); ncx = T1(nx - 1); ncy = T1(ny - 1); ncz = T1(nz - 1)
    nc = T1(ncx * ncy * ncz)
    _tmp_diff_vec_ = diff(x); all_equal = all(≈(_tmp_diff_vec_[1]), _tmp_diff_vec_)
    h1 = all_equal ? _tmp_diff_vec_[1] : error("grid spacing in x direction must be equal")
    _tmp_diff_vec_ = diff(y); all_equal = all(≈(_tmp_diff_vec_[1]), _tmp_diff_vec_)
    h2 = all_equal ? _tmp_diff_vec_[1] : error("grid spacing in y direction must be equal")
    _tmp_diff_vec_ = diff(z); all_equal = all(≈(_tmp_diff_vec_[1]), _tmp_diff_vec_)
    h3 = all_equal ? _tmp_diff_vec_[1] : error("grid spacing in z direction must be equal")
    h1 ≈ h2 ≈ h3 || error("Grid spacing in x, y and z directions must be equal")
    h    = T2(h1)
    invh = T2(1 / h)
    x1   = T2(x[1])
    x2   = T2(x[end])
    y1   = T2(y[1])
    y2   = T2(y[end])
    z1   = T2(z[1])
    z2   = T2(z[end])
    ζs   = T2(ζs)
    ms   = zeros(T2, ni)
    ps   = zeros(T2, ni, 3)
    vs   = zeros(T2, ni, 3)
    vsT  = zeros(T2, ni, 3)
    fs   = zeros(T2, ni, 3)
    vsxi = Array{T1}(vsxi)
    vsxv = Array{T2}(vsxv)
    vsyi = Array{T1}(vsyi)
    vsyv = Array{T2}(vsyv)
    vszi = Array{T1}(vszi)
    vszv = Array{T2}(vszv)
    ext  = ext == [0] ? (_tmp1=rand(T2, 3, 3), _tmp2=rand(T2, 3, 3)) : ext
    return Grid{T1, T2, Array{T1, 1}, AbstractArray{T2, 1}, AbstractArray{T2, 2}, NamedTuple}(
        ni, nx, ny, nz, nc, ncx, ncy, ncz, h, invh, x1, x2, y1, y2, z1, z2, ζs, ms, ps, vs, 
        vsT, fs, vsxi, vsxv, vsyi, vsyv, vszi, vszv, ext)
end