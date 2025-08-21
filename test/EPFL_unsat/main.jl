#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : 3d_druckerprager.jl                                                        |
|  Description: Case used to vaildate the functions                                        |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
+==========================================================================================#

using MaterialPointGenerator
using MaterialPointSolver
using MaterialPointVisualizer
using WGLMakie

include(joinpath(@__DIR__, "funcs.jl"))

init_h     = 0.04
init_ϵ     = :double
init_mtl   = :linearelastic
init_basis = :linear
init_NIC   = 8
init_FLIP  = 1.0
init_G     = -9.8
init_ρs    = 2e3
init_ρw    = 1e3
init_ν     = 0.4
init_Es    = 1.3e6
init_Gs    = init_Es / (2 * (1 + init_ν))
init_Ks    = init_Es / (3 * (1 - 2 * init_ν))
init_Kw    = 2e9
init_n     = 0.2975
init_k     = 4.41e-6
init_S     = 1.0
init_dev   = :cpu
init_T     = 30 * 1e-5
init_ΔT    = 1e-5
init_h5    = floor(Int, init_T / init_ΔT / 50)
init_var   = (:ξ,:Ω, :ρs, (:σw,))

# args setup
conf = init_conf(dev=init_dev, Δt=init_ΔT, t_tol=init_T,# h5_int=init_h5, varnames=init_var,
    prjpath=@__DIR__, prjname="EPFL_unsat", basis=init_basis, material=init_mtl)

# grid and boundary conditions setup
rangex  = -2*init_h:init_h:0.04+2*init_h
rangey  = -2*init_h:init_h:0.04+2*init_h
rangez  = -2*init_h:init_h:1.00+2*init_h
bg      = gridbuilder(rangex, rangey, rangez)
tmp_idx = findall(i -> bg.ξ[i, 1]≤0 || bg.ξ[i, 1]≥0.04, 1:bg.ni)
tmp_idy = findall(i -> bg.ξ[i, 2]≤0 || bg.ξ[i, 2]≥0.04, 1:bg.ni)
tmp_idz = findall(i -> bg.ξ[i, 3]≤0, 1:bg.ni)
tmp_idp = findall(i -> bg.ξ[i, 3]>1-init_h, 1:bg.ni)
vx_idx  = zeros(Int, bg.ni); vx_idx[tmp_idx] .= 1
vy_idx  = zeros(Int, bg.ni); vy_idx[tmp_idy] .= 1
vz_idx  = zeros(Int, bg.ni); vz_idx[tmp_idz] .= 1
vp_idx  = zeros(Int, bg.ni); vp_idx[tmp_idp] .= 1
grid    = init_grid(rangex, rangey, rangez, NIC=init_NIC,
    vsxi=vx_idx, vsxv=zeros(bg.ni), vsyi=vy_idx, vsyv=zeros(bg.ni), vszi=vz_idx, vszv=zeros(bg.ni),
    ext = (topS = [1.0], 
           mw   = zeros(Float64, bg.ni),
           mi   = zeros(Float64, bg.ni), 
           pw   = zeros(Float64, bg.ni, 3),
           fw   = zeros(Float64, bg.ni, 3),
           vw   = zeros(Float64, bg.ni, 3),
           vwT  = zeros(Float64, bg.ni, 3),
           vwxi = vx_idx, vwxv = zeros(bg.ni), 
           vwyi = vy_idx, vwyv = zeros(bg.ni), 
           vwzi = vp_idx, vwzv = zeros(bg.ni))
)

# material point setup
mh = grid.h * 0.5
ξ0 = meshbuilder(mh*0.5 : mh : 0.04-mh*0.5,
                 mh*0.5 : mh : 0.04-mh*0.5,
                 mh*0.5 : mh : 1.00-mh*0.5)
np = size(ξ0, 1)
nid = ones(Int, np)
mpts = init_mpts(ϵ=init_ϵ, NIC=init_NIC,
    ξ    = ξ0,
    h    = mh,
    G    = init_G,
    FLIP = init_FLIP,
    nid  = nid,
    ρs   = fill(init_ρs, np),
    Ks   = [init_Ks],
    Es   = [init_Es],
    Gs   = [init_Gs],
    ext = (
        σw = zeros(Float64, np),
        vw = zeros(Float64, np, 3),
        k  = init_k,
        S  = fill(init_S, np),
        kr = fill(kr(init_S), np),
        n  = fill(init_n, np),
        ρw = init_ρw,
        Kw = init_Kw,)
)
K0 = init_ν / (1 - init_ν)
γ0 = (1-init_n) * (init_ρs - init_S * init_ρw) * 9.8
@. mpts.σij[:, 1] = K0 * γ0 * (1-mpts.ξ[:, 3])
@. mpts.σij[:, 2] = K0 * γ0 * (1-mpts.ξ[:, 3])
@. mpts.σij[:, 1] =      γ0 * (1-mpts.ξ[:, 3])                               


# solver setup
vpos = Observable(mpts.ξ)
vcol = Observable(mpts.ext.σw)
let     
    custom_dark = theme_dark()
    custom_dark.textcolor = :white
    custom_dark.linecolor = :white
    custom_dark.backgroundcolor = "#181818"
    Makie.set_theme!(WGLMakie=(resize_to=:body,), custom_dark)
    fig = Figure()
    ax = LScene(fig[1, 1], show_axis=false)
    #canvas = ax.scene.plots[1]
    #canvas.ticks[:textcolor] = :white
    #canvas.frame[:axiscolor] = "#818181"
    #canvas.names[:textcolor] = :white
    p1 = scatter!(ax, vpos, color=vcol, colormap=:jet,
        inspector_label = (self, i, pos) -> (
        "x=" * string(pos[1]) *
        "\ny=" * string(pos[2]) *
        "\nz=" * string(pos[3]) *
        "\nval=" * string(vcol[][i])
    ))
    Colorbar(fig[1, 2], p1)
    DataInspector(fig, textcolor=:black,)

    display(fig)
end


mpmsolver!(tprocedure!, conf, grid, mpts, vpos, vcol)
#h5conf = (prjdst=conf.prjdst, prjname=conf.prjname)
#animation(h5conf)