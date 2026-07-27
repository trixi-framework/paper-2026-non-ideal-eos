const FIGDIR = joinpath(@__DIR__, "figs")
if !isdir(FIGDIR)
    mkpath(FIGDIR)
end

using Trixi
using Plots
using LaTeXStrings
include("trixi_utils.jl")

# run to define sol/semi before the loop (avoids world-age issues)
include("trixi_transcritical_shock.jl")  

for (polydeg, initial_refinement_level) in ((3, 7), (7, 6))
    # plot APEC-EC vs APEC-only solution
    trixi_include("trixi_transcritical_shock.jl",
                    volume_flux = flux_terashima_etal,
                    volume_integral = :(volume_integral_apec),
                    polydeg = polydeg,
                    initial_refinement_level = initial_refinement_level)

    # density
    pd = PlotData1D(sol.u[end], semi)
    p1 = plot(pd["rho"], linewidth=3, label="APEC only", dpi=400)

    # pressure
    p2 = plot(pd["p"],  linewidth=3, label="APEC only", dpi=400)

    trixi_include("trixi_transcritical_shock.jl",
                volume_flux = flux_terashima_etal,
                polydeg = polydeg,
                initial_refinement_level = initial_refinement_level)

    # u, x = sol_and_coordinates(sol)
    # xv = x[[1,end],:]
    # xc = vec(sum(x, dims = 1) / size(x, 1))
    # dx = xc[2] - xc[1]
    # w = semi.solver.basis.weights
    # u_avg = vec(w' * u / sum(w))
    # p_avg = vec(w' * pressure.(u, equations) / sum(w))

    # density
    pd = PlotData1D(sol.u[end], semi)
    plot!(p1, pd["rho"], linewidth=2, label="APEC + EC")
    # scatter!(p1, xc, getindex.(u_avg, 1), ms = 2, label="Cell average")
    plot!(p1, legendfontsize=12, guidefontsize=14, tickfontsize=12, titlefontsize=14, legend=:topright, dpi=400)
    savefig(p1, joinpath(FIGDIR, "transcritical_shock_rho_polydeg_$(polydeg).png"))

    # pressure
    plot!(p2, pd["p"], linewidth=2, label="APEC + EC")
    # scatter!(p2, xc, p_avg, ms = 2, label="Cell average")
    plot!(p2, legendfontsize=12, guidefontsize=14, tickfontsize=12, titlefontsize=14, legend=:topright, dpi=400)
    savefig(p2, joinpath(FIGDIR, "transcritical_shock_p_polydeg_$(polydeg).png"))
end

# blending coefficient plots
for volume_flux in (flux_central, flux_terashima_etal, flux_central_terashima_etal)
    trixi_include("trixi_transcritical_shock.jl",
                volume_flux = volume_flux,
                polydeg = 3,
                initial_refinement_level = 7)

    u, x = sol_and_coordinates(sol)
    xv = x[[1,end],:]
    xc = vec(sum(x, dims = 1) / size(x, 1))
    dx = xc[2] - xc[1]

    # density on left axis
    pd = PlotData1D(sol.u[end], semi)
    p3 = plot(pd["rho"], linewidth=2, label=L"$\rho$", dpi=400)

    # rerun RHS to refresh indicator alpha
    fill!(indicator_ec.cache.alpha, 0.0)
    du = similar(sol.u[end])
    Trixi.rhs!(du, sol.u[end], semi, sol.t[end])
    alpha = vec(indicator_ec.cache.alpha)

    bar!(twinx(), xc, alpha;
        label=nothing,
        alpha=0.45,          # bar transparency (not the indicator alpha)
        color=:orange,
        bar_width=1.1 * dx,  # one bar per cell
        bar_edges=false,
        legend=:topright,
        xlims=[-0.5, 0.5],
        ylims=[0.0, 0.0425])
    plot!(legendfontsize=12, guidefontsize=16, tickfontsize=14, dpi=400,
          titlefontsize=14, legend=false)

    savefig(p3, joinpath(FIGDIR, "transcritical_shock_indicator_$(volume_flux).png"))
end
# plot(p1, p2)
