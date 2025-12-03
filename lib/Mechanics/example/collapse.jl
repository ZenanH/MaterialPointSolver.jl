#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: 3D granular collapse test                                                  |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
+==========================================================================================#

using MaterialPointGenerator
using MaterialPointSolver
using MaterialPointVisualizer
using CUDA

using MaterialPointSolver.Mechanics # import the extension, `mpm_mechanics!` is from this extension

init_h     = 0.0025
init_ϵ     = :double
init_FLIP  = 1.0
init_ζs    = 0.0
init_G     = -9.8
init_ρs    = 2650
init_ν     = 0.3
init_Ks    = 7e5
init_Es    = init_Ks * (3 * (1 - 2 * init_ν))
init_Gs    = init_Es / (2 * (1 +     init_ν))
init_σt    = 0.0
init_ϕ     = deg2rad(19.8)
init_dev   = :cuda
init_T     = 0.6
init_Tcur  = 0.0
init_ΔT    = 0.5 * init_h / sqrt(init_Es / init_ρs)
init_h5    = 200
init_var   = (:ξ, :ϵq, :σij)

# args setup
conf = init_conf(dev=init_dev, Δt=init_ΔT, t_tol=init_T, h5_int=init_h5, varnames=init_var,
    prjpath=joinpath(@__DIR__, ".output"), prjname="Collapse")

# grid and boundary conditions setup
margin  = 0.
rangex  = [-0.02-margin, 0.07+margin]
rangey  = [-0.02-margin, 0.75+margin]
rangez  = [-0.02-margin, 0.12+margin]
bg      = gridbuilder(rangex, rangey, rangez, init_h)
tmp_idx = findall(i -> bg.ξ[i, 1]≤0||bg.ξ[i, 1]≥0.05||bg.ξ[i, 3]≤0||bg.ξ[i, 2]≤0, 1:bg.ni)
tmp_idy = findall(i -> bg.ξ[i, 2]≤0||bg.ξ[i, 3]≤0.00, 1:bg.ni)
tmp_idz = findall(i -> bg.ξ[i, 3]≤0, 1:bg.ni)
vx_idx  = zeros(Int, bg.ni); vx_idx[tmp_idx] .= 1
vy_idx  = zeros(Int, bg.ni); vy_idx[tmp_idy] .= 1
vz_idx  = zeros(Int, bg.ni); vz_idx[tmp_idz] .= 1 
grid    = init_grid(bg, ϵ=init_ϵ, ζs=init_ζs,
    vsxi=vx_idx, vsxv=zeros(bg.ni), 
    vsyi=vy_idx, vsyv=zeros(bg.ni), 
    vszi=vz_idx, vszv=zeros(bg.ni))

# material point setup
mh = grid.h * 0.5
ξ0 = meshbuilder([mh*0.5, 0.05-mh*0.5], 
                 [mh*0.5, 0.20-mh*0.5], 
                 [mh*0.5, 0.10-mh*0.5], mh) # xrange: 0:0.05, yrange: 0:0.2, zrange: 0:0.1
np = size(ξ0, 1)
mpts = init_mpts(ϵ=init_ϵ,
    ξ    = ξ0,
    h    = mh,
    G    = init_G,
    FLIP = init_FLIP,
    nid  = ones(Int, np),
    ρs   = fill(init_ρs, np),
    Ks   = [init_Ks],
    Es   = [init_Es],
    Gs   = [init_Gs],
    ϕ    = [init_ϕ]
)

# solver setup
mpmsolver!(mpm_mechanics!, conf, grid, mpts)
h5conf = (prjdst=conf.prjdst, prjname=conf.prjname)
animation(h5conf)