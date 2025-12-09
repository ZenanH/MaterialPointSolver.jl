export Particle, DeviceParticle
export init_mpts

abstract type DeviceParticle{T1, T2} end

struct Particle{T1, T2,
    T3 <: AbstractArray{T1, 1},
    T4 <: AbstractArray{T2, 1},
    T5 <: AbstractArray{T2, 2},
    T6 <: NamedTuple
} <: DeviceParticle{T1, T2}
    np  ::T1
    h   ::T2
    G   ::T2
    FLIP::T2
    nid ::T3
    cfl ::T4
    ξ   ::T5
    ξ0  ::T5
    F   ::T5
    ext ::T6
end

@KAadapt Particle

function init_mpts(; 
    ϵ::Symbol=:double, 
    ξ::AbstractArray, 
    h::Real, G::Real=-9.8, FLIP::Real=1.0, 
    nid::AbstractArray=[0], 
    ext::NamedTuple=(_a=0, _b=0)
)
    T1   = ϵ == :single ? Int32   : Int64
    T2   = ϵ == :single ? Float32 : Float64
    np   = T1(size(ξ, 1))
    h    = T2(h)
    G    = T2(G)
    FLIP = 0 ≤ T2(FLIP) ≤ 1 ? T2(FLIP) : error("FLIP must be in [0, 1]")
    nid  = nid == [0] ? ones(T1, np) : nid
    cfl  = zeros(T2, np)
    ξ    = Array{T2}(ξ)
    ξ0   = copy(ξ)
    F    = Array{T2}(repeat([1 0 0 0 1 0 0 0 1], np))
    ext  = ext == (_a=0, _b=0) ? (_tmp1=rand(T2, 3, 3), _tmp2=rand(T2, 3, 3)) : ext
    return KAupload(Array, Particle{T1, T2, Array{T1, 1}, Array{T2, 1}, Array{T2, 2}, NamedTuple}(
        np, h, G, FLIP, nid, cfl, ξ, ξ0, F, ext))
end