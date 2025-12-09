export Grid, DeviceGrid
export init_grid

abstract type DeviceGrid{T1, T2} end

struct Grid{T1, T2, T3 <: NamedTuple} <: DeviceGrid{T1, T2}
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
    ext ::T3
end

@KAadapt Grid

function init_grid(bg::NamedTuple; ϵ::Symbol=:double, ζs=0.0, ext=(_a=0, _b=0))
    # Determine types
    T1  = ϵ == :single ? Int32   : Int64
    T2  = ϵ == :single ? Float32 : Float64
    # Initialize fields
    nx = T1(bg.nx); ncx = T1(nx - 1)
    ny = T1(bg.ny); ncy = T1(ny - 1)
    nz = T1(bg.nz); ncz = T1(nz - 1)
    ni = T1(bg.ni); nc = T1(ncx * ncy * ncz)
    nc == bg.nc || error("Inconsistent grid size!")
    h    = T2(bg.h)
    invh = T2(bg.inv_h)
    x1   = T2(bg.x1)
    x2   = T2(bg.x2)
    y1   = T2(bg.y1)
    y2   = T2(bg.y2)
    z1   = T2(bg.z1)
    z2   = T2(bg.z2)
    ζs   = T2(ζs)
    ext  = ext == (_a=0, _b=0) ? (_tmp1=rand(T2, 3, 3), _tmp2=rand(T2, 3, 3)) : ext
    return KAupload(Array, Grid{T1, T2, NamedTuple}(ni, nx, ny, nz, nc, ncx, ncy, ncz, h, 
        invh, x1, x2, y1, y2, z1, z2, ζs, ext))
end