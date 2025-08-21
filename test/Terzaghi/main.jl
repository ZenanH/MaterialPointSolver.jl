#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: Terzaghi consolidation test (saturated MPM)                                |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
+==========================================================================================#

using MaterialPointGenerator
using MaterialPointSolver
using MaterialPointVisualizer
# using WGLMakie

include(joinpath(@__DIR__, "funcs.jl"))

init_h     = 0.05
init_ϵ     = :double
init_mtl   = :linearelastic
init_basis = :linear
init_NIC   = 8
init_FLIP  = 1.0
init_G     = 0.0
init_ρs    = 2.65e3
init_ρw    = 1e3
init_ν     = 0.0
init_Es    = 1e7
init_Gs    = init_Es / (2 * (1 + init_ν))
init_Ks    = init_Es / (3 * (1 - 2 * init_ν))
init_Kw    = 2.2e9
init_n     = 0.4
init_k     = 1e-3
init_S     = 1.0
init_σw    = -1e4
init_dev   = :cpu
init_T     = 1.0
init_Tcur  = 0.0
init_ΔT    = 1e-5
init_h5    = floor(Int, init_T / init_ΔT / 200)
init_var   = (:ξ, (:σw,))

# args setup
conf = init_conf(dev=init_dev, Δt=init_ΔT, t_tol=init_T, h5_int=init_h5, varnames=init_var,
    prjpath=@__DIR__, prjname="Terzaghi", basis=init_basis, material=init_mtl)

# grid and boundary conditions setup
rangex  = -2*init_h:init_h:0.1+2*init_h
rangey  = -2*init_h:init_h:0.1+2*init_h
rangez  = -2*init_h:init_h:1.0+2*init_h
bg      = meshbuilder(rangex, rangey, rangez)
ni      = size(bg, 1)
tmp_idx = findall(i -> bg[i, 1]≤0||bg[i, 1]≥0.1, 1:ni)
tmp_idy = findall(i -> bg[i, 2]≤0||bg[i, 2]≥0.1, 1:ni)
tmp_idz = findall(i -> bg[i, 3]≤0, 1:ni)
vx_idx  = zeros(Int, ni); vx_idx[tmp_idx] .= 1
vy_idx  = zeros(Int, ni); vy_idx[tmp_idy] .= 1
vz_idx  = zeros(Int, ni); vz_idx[tmp_idz] .= 1 
grid    = init_grid(rangex, rangey, rangez, NIC=init_NIC,
    vsxi=vx_idx, vsxv=zeros(ni), vsyi=vy_idx, vsyv=zeros(ni), vszi=vz_idx, vszv=zeros(ni),
    ext = (mw   = zeros(Float64, ni),
           mi   = zeros(Float64, ni),
           pw   = zeros(Float64, ni, 3),
           fw   = zeros(Float64, ni, 3),
           vw   = zeros(Float64, ni, 3),
           vwT  = zeros(Float64, ni, 3),
           vwxi = vx_idx,
           vwxv = zeros(ni),
           vwyi = vy_idx,
           vwyv = zeros(ni),
           vwzi = vz_idx,
           vwzv = zeros(ni))
)

# material point setup
mh = grid.h * 0.5
ξ0 = meshbuilder(mh*0.5 : mh : 0.1-mh*0.5,
                 mh*0.5 : mh : 0.1-mh*0.5,
                 mh*0.5 : mh : 1.0-mh*0.5)
np = size(ξ0, 1)
bcp = findall(i -> ξ0[i, 3] ≥ 1 - mh, 1:np)
mpts = init_mpts(ϵ=init_ϵ, NIC=init_NIC,
    ξ    = ξ0,
    h    = mh,
    G    = init_G,
    FLIP = init_FLIP,
    nid  = ones(Int, np),
    ρs   = fill(init_ρs, np),
    Ks   = [init_Ks],
    Es   = [init_Es],
    Gs   = [init_Gs],
    ext = (
        σw  = fill(init_σw, np),
        vw  = zeros(Float64, np, 3),
        k   = init_k,
        kr  = fill(1.0, np),
        S   = fill(init_S, np),
        n   = fill(init_n, np),
        Kw  = init_Kw,
        ρw  = init_ρw,
        bcp = bcp,
        mpf = -100 / length(bcp))
)

# solver setup
mpmsolver!(tprocedure!, conf, grid, mpts)
h5conf = (prjdst=conf.prjdst, prjname=conf.prjname)
animation(h5conf)

# let
#     custom_dark = theme_dark()
#     custom_dark.textcolor = :white
#     custom_dark.linecolor = :white
#     custom_dark.backgroundcolor = "#181818"
#     Makie.set_theme!(WGLMakie=(resize_to=:body,), custom_dark)
#     fig = Figure()
#     ax = LScene(fig[1, 1], show_axis=true)
#     canvas = ax.scene.plots[1]
#     canvas.ticks[:textcolor] = :white
#     canvas.ticks[:fontsize] = (15, 15, 15)
#     canvas.frame[:axiscolor] = "#818181"
#     canvas.names[:textcolor] = :white
#     canvas.names[:fontsize] = (15, 15, 15)
#     scatter!(ax, mpts.ξ, color=mpts.ext.σw, colormap=:turbo,
#         inspector_label = (self, i, pos) -> (
#         "x=" * string(pos[1]) *
#         "\ny=" * string(pos[2]) *
#         "\nz=" * string(pos[3]) *
#         "\nval=" * string(mpts.ext.σw[i])
#     ))
#     scatter!(ax, bg, markersize=4, alpha=0.2)
#     #scatter!(ax, mpts.ξ[bcp, :])
#     #scatter!(ax, bg[findall(grid.ext.vwyi.==1), :],)
#     DataInspector(fig, textcolor = :black,)
#     display(fig)
# end