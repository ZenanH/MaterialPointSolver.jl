using CairoMakie
using MaterialPointSolver
using Symbolics


@variables Pc # define the variable name
S = 1 - 0.10152 * (Pc/9.8e3)^2.4279 # define the equation
dSdPc = Symbolics.derivative(S, Pc) # calculate the derivative

let
    figfontr = MaterialPointSolver.tnr
    figfontb = MaterialPointSolver.tnrb
    figfonti = MaterialPointSolver.tnri
    fig = Figure(size=(400, 270), fonts=(; regular=figfontr, bold=figfontb, italic=figfonti))
    ax1 = Axis(fig[1, 1], ylabel="Degree of saturation", xlabel="Matric suction (Pa)")
    ax2 = Axis(fig[1, 1], ylabel="Hydraulic conductivity (mm/s)")
    ax2.yaxisposition = :right
    ax2.flip_ylabel = true
    hidespines!(ax2)
    hidexdecorations!(ax2)
    P = range(0, 9200, length=200)
    swcc(s) = 1 - 0.10152 * (s/9.8e3)^2.4279
    hcc(s) = 1.0 - 2.207 * (1 - swcc(s))
    SWCC_values = [swcc(s) for s in P]
    HCC_values = [hcc(s) for s in P] .* 1e3 .* 4.41e-6
    p1 = lines!(ax1, hcat(P, SWCC_values), color=:blue)
    p2 = lines!(ax2, hcat(P, HCC_values), color=:red, linestyle=:dashdot, linewidth=2)
    axislegend(ax1, [p1, p2], ["SWCC", "HCC"], position=:lb)
    display(fig)
end

let
    figfontr = MaterialPointSolver.tnr
    figfontb = MaterialPointSolver.tnrb
    figfonti = MaterialPointSolver.tnri
    fig = Figure(size=(400, 270), fonts=(; regular=figfontr, bold=figfontb, italic=figfonti))
    ax1 = Axis(fig[1, 1], ylabel="∂S/∂P", xlabel="Matric suction (Pa)")
    P = range(0, 9200, length=200)
    ∂swcc(s) = -5.03e-11(s^1.4279)
    ∂SWCC_values = [∂swcc.(p) for p in P]
    p1 = lines!(ax1, hcat(P, ∂SWCC_values), color=:blue)
    display(fig)
end