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

function init_grid(bg::NamedTuple; 
    ϵ::Symbol=:double, vsxi=[1], vsxv=[0], vsyi=[1], vsyv=[0], vszi=[1], vszv=[0], 
    ζs=0.0, ext=[0]
)
    T1  = ϵ == :single ? Int32   : Int64
    T2  = ϵ == :single ? Float32 : Float64

    nx = T1(bg.nx); ncx = T1(nx - 1)
    ny = T1(bg.ny); ncy = T1(ny - 1)
    nz = T1(bg.nz); ncz = T1(nz - 1)
    ni = T1(bg.ni); nc = T1(ncx * ncy * ncz); nc == bg.nc
    h    = T2(bg.h)
    invh = T2(bg.inv_h)
    x1   = T2(bg.x1)
    x2   = T2(bg.x2)
    y1   = T2(bg.y1)
    y2   = T2(bg.y2)
    z1   = T2(bg.z1)
    z2   = T2(bg.z2)
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