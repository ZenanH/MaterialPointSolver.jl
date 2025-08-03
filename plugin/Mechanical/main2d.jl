using OpenMPM
using OpenMPM.Preprocessing
using OpenMPM.Postprocessing
using Metal
#using CUDA
include(joinpath(@__DIR__, "mechanical.jl"))

init_h     = 0.0025
init_ϵ     = set_precision(:single)
init_basis = Linear()
init_dim   = Dim2D()
init_FLIP  = 1.0
init_G     = -9.8
init_ρs    = 2650
init_ν     = 0.3
init_Ks    = 7e5
init_Es    = init_Ks * (3 * (1 - 2 * init_ν))
init_Gs    = init_Es / (2 * (1 +     init_ν))
init_σt    = 0.0
init_ϕ     = deg2rad(19.8)
init_ψ     = 0.0
init_c     = 0.0
init_dev   = :metal
init_T     = 1
init_Tcur  = 0.0
init_ΔT    = 0.5 * init_h / sqrt(init_Es / init_ρs)
init_h5    = hasHDF5(interval=floor(Int, init_T / init_ΔT / 200), varnames=(:ξ, :ϵq))

# args setup
conf = generate_config(init_ϵ, init_dim, init_basis, init_dev, init_h5, t_tol=init_T, 
    Δt=init_ΔT, prjpath=@__DIR__, prjname="2D_example")

# grid and boundary conditions setup
bg      = gridbuilder(-0.025:init_h:0.82, -0.025:init_h:0.12)
tmp_idx = findall(i -> bg.ξ[1, i]≤0.0||bg.ξ[1, i]≥0.8||bg.ξ[2, i]≤0, 1:bg.ni)
tmp_idy = findall(i -> bg.ξ[2, i]≤0.0, 1:bg.ni)
vx_idx  = zeros(Int, bg.ni); vx_idx[tmp_idx] .= 1
vy_idx  = zeros(Int, bg.ni); vy_idx[tmp_idy] .= 1 
grid    = generate_fields(init_ϵ, -0.025:init_h:0.82, -0.025:init_h:0.12, 
    ms       = zeros(bg.ni),
    Ω        = zeros(bg.ni),
    vs       = zeros(2, bg.ni),
    ps       = zeros(2, bg.ni),
    fs       = zeros(2, bg.ni),
    vsT      = zeros(2, bg.ni),
    vx_s_idx = vx_idx,
    vx_s_val = zeros(bg.ni),
    vy_s_idx = vy_idx,
    vy_s_val = zeros(bg.ni)
)

# material point setup
mh = grid.h * 0.5
ξ0 = meshbuilder(mh*0.5 : mh : 0.2-mh*0.5, mh*0.5 : mh : 0.1-mh*0.5) # xrange: 0:0.2, yrange: 0:0.1
np = size(ξ0, 2)
mpts = generate_fields(init_ϵ,
    ξ0   = ξ0,
    ξ    = copy(ξ0),
    ρs   = ones(np) .* init_ρs,
    ρs0  = ones(np) .* init_ρs,
    Ω    = ones(np) .* mh^2,
    Ω0   = ones(np) .* mh^2,
    vs   = zeros(2, np),
    σij  = zeros(4, np),
    F    = repeat(float.([1, 0, 0, 1]), 1, np),
    ϵq   = zeros(np),
    ϵk   = zeros(np),
    np   = np,
    p2n  = zeros(Int, conf.NIC, np),
    Nij  = zeros(conf.NIC, np),
    ∂Nx  = zeros(conf.NIC, np),
    ∂Ny  = zeros(conf.NIC, np),
    nid  = ones(Int, np),
    FLIP = init_FLIP,
    G    = init_G,
    ν    = [init_ν],
    Es   = [init_Es],
    Gs   = [init_Gs],
    Ks   = [init_Ks],
    σt   = [init_σt],
    ϕ    = [init_ϕ],
    ψ    = [init_ψ],
    c    = [init_c],
    NIC  = conf.NIC,
)

mpmsolver!(Mechanical.procedure!, conf, grid, mpts)

animation(conf)

rm(conf.prjdst, force=true, recursive=true) # clean up project directory