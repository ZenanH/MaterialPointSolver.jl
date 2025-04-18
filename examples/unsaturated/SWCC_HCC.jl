using MaterialPointSolver
using CairoMakie

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
    S_min = 0.125
    S_max = 1.0
    λ = 0.7
    P_ref = 3e3
    k_sat = 4.71e-3
    α = 2e-4
    β = 2.214
    P = 10 .^ range(log10(1e2), log10(1e5), length=200)
    SWCC_values = [SWCC(S_min, S_max, p, P_ref, λ) for p in P] .* 100
    HCC_values = [HCC(k_sat, α, β, p) for p in P] .* 1e3
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
    S_min = 0.125
    S_max = 1.0
    λ     = 0.7
    P_ref = 3e3
    P     = 10 .^ range(log10(1e2), log10(1e5), length=200)
    ∂SWCC_values = ∂SWCC.(S_min, S_max, P, P_ref, λ)
    p1 = lines!(ax1, hcat(P, ∂SWCC_values), color=:blue)
    display(fig)
end