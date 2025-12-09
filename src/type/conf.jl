export Config
export HDF5Config, H5_T, H5_F
export dev_backend, init_conf

abstract type HDF5Config end

struct H5_T <: HDF5Config
    gname    ::Ref{Int}
    k        ::Ref{Int}
    tol_iters::Int
    interval ::Vector{Float64}
    varnames ::Tuple
    fpvar
end

struct H5_F <: HDF5Config
    iters::Ref{Int}
end

struct Config
    dev
    adaptive::Bool
    h5      ::HDF5Config
    iters   ::Ref{Int64}
    log_int ::Float64
    prjname ::String
    prjpath ::String
    prjdst  ::String
    stime   ::Ref{Float64}
    etime   ::Ref{Float64}
    Δt      ::Float64
    t_tol   ::Float64
    t_cur   ::Float64
    t_eld   ::Float64
    αT      ::Float64
end

dev_backend(sym::Symbol=:cpu) = dev_backend(Val(sym))
dev_backend(::Val{S}) where {S} = _missing_backend(S)
_missing_backend(S) = throw(ArgumentError(
    """
    Backend :$S is unavailable.
    Please 1) \e[1;36madd\e[0m and ensure compatibility with the corresponding packages:
    ────────────┬────────────┬──────────────
     GPU vendor │ Package    │ Backend name
    ────────────┼────────────┼──────────────
     Nvidia     │ \e[1;32mCUDA.jl\e[0m    │ \e[0;31m:cuda\e[0m
     AMD        │ \e[1;32mAMDGPU.jl\e[0m  │ \e[0;31m:rocm\e[0m
     Intel      │ \e[1;32moneAPI.jl\e[0m  │ \e[0;31m:oneapi\e[0m
     Apple      │ \e[1;32mMetal.jl\e[0m   │ \e[0;31m:metal\e[0m
     CPU        │ by default │ \e[0;31m:cpu\e[0m
    ────────────┴────────────┴──────────────
    and 2) \e[1;36mconfirm\e[0m that the Julia version \e[1;31m>= 1.9\e[0m to enable extensions.
    """
))
dev_backend(::Val{:cpu}) = CPU()

function init_conf(; dev::Symbol=:cpu, adaptive::Bool=false, h5_int::Int=0, 
    varnames::Tuple=(:default,), log_int::Real=3.0, prjname, prjpath, Δt, t_tol, 
    t_cur::Real=0.0, t_eld::Real=0.0, αT::Real=0.5
)
    dev = dev_backend(dev)
    # 展平字段路径：(:x, :y, (:u, :z)) → [(:x,), (:y,), (:ext, :u), (:ext, :z)]
    field_paths = Tuple[]
    for item in varnames
        if item isa Symbol
            push!(field_paths, (item,))
        elseif item isa Tuple
            for subfield in item
                push!(field_paths, (:ext, subfield))
            end
        else
            error("Unsupported entry in varnames: $item")
        end
    end

    h5 = h5_int > 0 ?
        H5_T(Ref(1), Ref(1), h5_int, collect(range(t_cur, t_tol; length=h5_int)), varnames, field_paths) :
        H5_F(Ref(0))

    iters = Ref{Int64}(0)
    prjdst = joinpath(prjpath, prjname)
    isdir(prjdst) && rm(prjdst; recursive=true, force=true); mkpath(prjdst)
    stime = Ref{Float64}(time())
    etime = Ref{Float64}(time())
    Δt    = Float64(Δt)
    t_tol = Float64(t_tol)
    t_cur = Float64(t_cur)
    t_eld = Float64(t_eld)
    αT    = Float64(αT)

    return Config(dev, adaptive, h5, iters, log_int, prjname, prjpath, prjdst, stime, etime, 
        Δt, t_tol, t_cur, t_eld, αT)
end
