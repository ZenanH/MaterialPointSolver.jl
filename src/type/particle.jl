export Particle, DeviceParticle
export init_mpts

abstract type DeviceParticle{T1, T2} end

struct Particle{T1, T2,
    T3 <: AbstractArray{T1, 1},
    T4 <: AbstractArray{T1, 2},
    T5 <: AbstractArray{T2, 1},
    T6 <: AbstractArray{T2, 2},
    T7 <: NamedTuple
} <: DeviceParticle{T1, T2}
    np  ::T1
    NIC ::T1
    h   ::T2
    G   ::T2
    FLIP::T2
    nid ::T3
    p2n ::T4
    ϵq  ::T5
    ϵk  ::T5
    Ω   ::T5
    Ω0  ::T5
    ρs  ::T5
    ρs0 ::T5
    cfl ::T5
    Ks  ::T5
    Es  ::T5
    Gs  ::T5
    c   ::T5
    cr  ::T5
    ϕ   ::T5
    ϕr  ::T5
    ψ   ::T5
    Hp  ::T5
    σt  ::T5
    ξ   ::T6
    ξ0  ::T6
    σij ::T6
    vs  ::T6
    Nij ::T6
    ∂Nx ::T6
    ∂Ny ::T6
    ∂Nz ::T6
    F   ::T6
    ext ::T7
end

@KAadapt Particle

function init_mpts(; ϵ::Symbol=:double, ξ, NIC=27, h, G, FLIP, nid=[0], ρs, Ks=[0.0], 
    Es=[0.0], Gs=[0.0], c=[0.0], cr=[0.0], ϕ=[0.0], ϕr=[0.0], ψ=[0.0], Hp=[0.0], σt=[0.0],
    ext=[0])
    T1   = ϵ == :single ? Int32   : Int64
    T2   = ϵ == :single ? Float32 : Float64
    np   = T1(size(ξ, 1))
    NIC  = T1(NIC)
    h    = T2(h)
    G    = T2(G)
    FLIP = T2(FLIP)
    nid  = nid == [0] ? ones(T1, np) : nid
    p2n  = zeros(T1, np, NIC) 
    ϵq   = zeros(T2, np)
    ϵk   = zeros(T2, np)
    Ω    = Array{T2}(fill(h^3, np))
    Ω0   = copy(Ω)
    ρs   = Array{T2}(ρs)
    ρs0  = copy(ρs)
    cfl  = zeros(T2, np)
    Ks   = Array{T2}(Ks)
    Es   = Array{T2}(Es)
    Gs   = Array{T2}(Gs)
    c    = Array{T2}(c)
    cr   = Array{T2}(cr)
    ϕ    = Array{T2}(ϕ)
    ϕr   = Array{T2}(ϕr)
    ψ    = Array{T2}(ψ)
    Hp   = Array{T2}(Hp)
    σt   = Array{T2}(σt)
    ξ    = Array{T2}(ξ)
    ξ0   = copy(ξ)
    σij  = zeros(T2, np, 6)
    vs   = zeros(T2, np, 3)
    Nij  = zeros(T2, np, NIC)
    ∂Nx  = zeros(T2, np, NIC)
    ∂Ny  = zeros(T2, np, NIC)
    ∂Nz  = zeros(T2, np, NIC)
    F    = Array{T2}(repeat([1 0 0 0 1 0 0 0 1], np))
    ext  = ext == [0] ? (_tmp1=rand(T2, 3, 3), _tmp2=rand(T2, 3, 3)) : ext
    return Particle{T1, T2, Array{T1, 1}, Array{T1, 2}, Array{T2, 1}, Array{T2, 2}, NamedTuple}(
        np, NIC, h, G, FLIP, nid, p2n, ϵq, ϵk, Ω, Ω0, ρs, ρs0, cfl, Ks, Es, Gs, c, cr, ϕ, 
        ϕr, ψ, Hp, σt, ξ, ξ0, σij, vs, Nij, ∂Nx, ∂Ny, ∂Nz, F, ext)
end