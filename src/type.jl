#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: Type system                                                                |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, ISTE, Université de Lausanne                                   |
|  Maintainer : Zenan Huo                                                                  |
+==========================================================================================#

export Precision, MPMDim, Dim2D, Dim3D
export MPMBasis, Linear, uGIMP, Bspline, Cspline
export MPMHDF5, hasHDF5, noHDF5
export set_precision, generate_fields, generate_config

struct Precision{T1, T2} end

abstract type MPMDim end
struct Dim2D <: MPMDim end
struct Dim3D <: MPMDim end

abstract type MPMBasis end
struct Linear  <: MPMBasis end
struct uGIMP   <: MPMBasis end
struct Bspline <: MPMBasis end
struct Cspline <: MPMBasis end

abstract type MPMHDF5 end
struct hasHDF5 <: MPMHDF5
    iters::Ref{Int}
    gname::Ref{Int}
    interval::Int
    varnames::Tuple{Vararg{Symbol}}
end
struct noHDF5 <:MPMHDF5 
    iters
end

function hasHDF5(; interval::Int, varnames::Tuple{Vararg{Symbol}})
    interval > 0 || throw(ArgumentError("interval must > 0."))
    return hasHDF5(Ref(0), Ref(1), interval, varnames)
end

noHDF5() = noHDF5(Ref{Int}(0))

function set_precision(ϵ::Symbol=:double)
    T1 = ϵ === :single ? Int32   : Int64
    T2 = ϵ === :single ? Float32 : Float64
    return Precision{T1, T2}()
end

@inline function _register_fields(ϵ::Precision{I, F}; kwargs...) where {I, F}
    isempty(kwargs) && throw(ArgumentError("At least one data field must be provided."))
    nt_in = NamedTuple(kwargs)
    new_vals = (
        x isa Signed        ? convert(I, x) :
        x isa AbstractFloat ? convert(F, x) :
        x isa AbstractArray ? begin
            T = eltype(x)
            if T <: Signed
                T == I ? x : convert.(I, x) # fast-path when already I
            elseif  T <: AbstractFloat
                T == F ? x : convert.(F, x) # fast-path when already F
            else
                throw(ArgumentError("Array element type $T is not <: AbstractFloat or Signed"))
            end
        end :
        throw(ArgumentError("Field type $(typeof(x)) is not supported"))
        for x in values(nt_in)
    ); return NamedTuple{keys(nt_in)}(new_vals)
end

generate_fields(ϵ::Precision; kwargs...) = _register_fields(ϵ; kwargs...)

function generate_fields(ϵ::Precision{I, F}, xr::AbstractRange, yr::AbstractRange; kwargs...) where {I, F}
    DISALLOWED_KEYS = (:nx, :ny, :ni, :nc, :h, :x1, :x2, :y1, :y2, :inv_h)
    @inbounds for k in keys(kwargs)
        k in DISALLOWED_KEYS && throw(ArgumentError("keyword $k is not allowed."))
    end; part1 = _register_fields(ϵ; kwargs...)
    x, y = convert.(Float64, xr), convert.(Float64, yr)
    nx = length(x); nx ≤ 1 && error("Grid must have at least 2 points in x-direction")
    ny = length(y); ny ≤ 1 && error("Grid must have at least 2 points in y-direction")
    ni = nx * ny; nc = (nx - 1) * (ny - 1)
    _tmp_diff_vec_ = diff(x); all_equal = all(≈(_tmp_diff_vec_[1]), _tmp_diff_vec_)
    h1 = all_equal ? _tmp_diff_vec_[1] : error("grid spacing in x direction must be equal")
    _tmp_diff_vec_ = diff(y); all_equal = all(≈(_tmp_diff_vec_[1]), _tmp_diff_vec_)
    h2 = all_equal ? _tmp_diff_vec_[1] : error("grid spacing in y direction must be equal")
    h1 ≈ h2 || error("Grid spacing in x and y directions must be equal")
    h = h1; inv_h = 1 / h
    part2 = _register_fields(ϵ; nx=nx, ny=ny, ni=ni, nc=nc, h=h, x1=x[1], x2=x[end], y1=y[1], y2=y[end], inv_h=inv_h)
    return merge(part1, part2)
end

function generate_fields(ϵ::Precision{I, F}, xr::AbstractRange, yr::AbstractRange, zr::AbstractRange; kwargs...) where {I, F}
    DISALLOWED_KEYS = (:nx, :ny, :nz, :ni, :nc, :h, :x1, :x2, :y1, :y2, :z1, :z2, :inv_h)
    @inbounds for k in keys(kwargs)
        k in DISALLOWED_KEYS && throw(ArgumentError("keyword $k is not allowed."))
    end; part1 = _register_fields(ϵ; kwargs...)
    x, y, z = convert.(Float64, xr), convert.(Float64, yr), convert.(Float64, zr)
    nx = length(x); nx ≤ 1 && error("Grid must have at least 2 points in x-direction")
    ny = length(y); ny ≤ 1 && error("Grid must have at least 2 points in y-direction")
    nz = length(z); nz ≤ 1 && error("Grid must have at least 2 points in z-direction")
    ni = nx * ny * nz; nc = (nx - 1) * (ny - 1) * (nz - 1)
    _tmp_diff_vec_ = diff(x); all_equal = all(≈(_tmp_diff_vec_[1]), _tmp_diff_vec_)
    h1 = all_equal ? _tmp_diff_vec_[1] : error("grid spacing in x direction must be equal")
    _tmp_diff_vec_ = diff(y); all_equal = all(≈(_tmp_diff_vec_[1]), _tmp_diff_vec_)
    h2 = all_equal ? _tmp_diff_vec_[1] : error("grid spacing in y direction must be equal")
    _tmp_diff_vec_ = diff(z); all_equal = all(≈(_tmp_diff_vec_[1]), _tmp_diff_vec_)
    h3 = all_equal ? _tmp_diff_vec_[1] : error("grid spacing in z direction must be equal")
    h1 ≈ h2 ≈ h3 || error("Grid spacing in x, y and z directions must be equal")
    h = h1; inv_h = 1 / h
    part2 = _register_fields(ϵ; nx=nx, ny=ny, nz=nz, ni=ni, nc=nc, h=h, x1=x[1], x2=x[end], y1=y[1], y2=y[end], z1=z[1], z2=z[end], inv_h=inv_h)
    return merge(part1, part2)
end

function generate_config(ϵ::Precision{I, F}, dim::MPMDim, basis::MPMBasis, dev::Symbol, 
    h5::MPMHDF5; t_tol::Real, t_cur::Real=0.0, Δt::Real, t_load::Real=0, 
    log_interval::Real=3.0, prjpath::String, prjname::String, kwargs...
) where {I, F}
    # input validation
    prjdst = joinpath(prjpath, prjname)
    isdir(prjdst) && rm(prjdst; recursive=true, force=true); mkpath(prjdst)
    dev in [:cpu, :cuda, :metal, :rocm, :oneapi] || throw(ArgumentError(
        "Device $dev is not supported. Supported devices are: :cpu, :cuda, :metal, :rocm, :oneapi."))
    t_tol > 0 || throw(ArgumentError("t_tol must be positive."))
    t_cur < t_tol || throw(ArgumentError("t_cur must be less than t_tol."))
    Δt > 0 || throw(ArgumentError("Δt must be positive."))
    0 ≤ t_load ≤ t_tol || throw(ArgumentError("t_load must be in [0, t_tol]."))
    log_interval > 0 || throw(ArgumentError("log_interval must be positive."))
    # convert types
    t_tol         = convert(F, t_tol)
    t_cur         = convert(F, t_cur)
    Δt            = convert(F, Δt)
    t_load        = convert(F, t_load)
    # new fields
    iters      = Ref{I}(0)
    time_start = Ref{F}(time())
    timer1     = Ref{F}(time())
    timer2     = Ref{F}(time())
    time_end   = Ref{F}(time())

    dev = dev_backend(dev)
    NIC = if dim isa Dim2D
        if basis isa Linear
            4
        elseif basis isa uGIMP || basis isa Bspline
            9
        elseif basis isa Cspline
            16
        end
    elseif dim isa Dim3D
        if basis isa Linear
            8
        elseif basis isa uGIMP || basis isa Bspline
            27
        elseif basis isa Cspline
            64
        end
    end

    return (; ϵ=ϵ, dim=dim, basis=basis, dev=dev, h5=h5, iters=iters, time_start=time_start, 
        timer1=timer1, timer2=timer2, time_end=time_end, t_tol=t_tol, t_cur=t_cur, Δt=Δt, 
        t_load=t_load, log_interval=log_interval, NIC=NIC, prjdst=prjdst, prjpath=prjpath, 
        prjname=prjname, kwargs...)
end