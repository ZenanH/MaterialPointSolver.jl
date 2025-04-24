using CairoMakie
using MaterialPointSolver

let
    figfontr = MaterialPointSolver.tnr
    figfontb = MaterialPointSolver.tnrb
    figfonti = MaterialPointSolver.tnri
    fig = Figure(size=(400, 270), fonts=(; regular=figfontr, bold=figfontb, italic=figfonti))
    ax1 = Axis(fig[1, 1], xscale=log10, ylabel="Degree of saturation (%)", xlabel="Matric suction (Pa)")
    ax2 = Axis(fig[1, 1], xscale=log10, ylabel="Hydraulic conductivity (mm/s)", xlabel="Suction (Pa)")
    ax2.yaxisposition = :right
    ax2.flip_ylabel = true
    hidespines!(ax2)
    hidexdecorations!(ax2)
    P = 10 .^ range(log10(2000), log10(18000), length=200)
    swcc(s) = 1 - 0.10152 * (s*0.000102)^2.4279
    hcc(s) = 1.0 - 2.207 * (1 - swcc(s))
    SWCC_values = [swcc(s) for s in P] .* 100
    HCC_values = [hcc(s) for s in P] .* 1e3 .* 4.41e-6
    p1 = lines!(ax1, hcat(P, SWCC_values), color=:blue)
    p2 = lines!(ax2, hcat(P, HCC_values), color=:red, linestyle=:dashdot, linewidth=2)
    axislegend(ax1, [p1, p2], ["SWCC", "HCC"], position=:rt)
    display(fig)
end

let
    figfontr = MaterialPointSolver.tnr
    figfontb = MaterialPointSolver.tnrb
    figfonti = MaterialPointSolver.tnri
    fig = Figure(size=(400, 270), fonts=(; regular=figfontr, bold=figfontb, italic=figfonti))
    ax1 = Axis(fig[1, 1], xscale=log10, ylabel="∂S/∂P", xlabel="Matric suction (Pa)")
    P = 10 .^ range(log10(2000), log10(18000), length=200)
    ∂swcc(s) = -5.0241518982722e-11*s^1.4279
    ∂SWCC_values = [∂swcc.(p) for p in P]
    p1 = lines!(ax1, hcat(P, ∂SWCC_values), color=:blue)
    display(fig)
end