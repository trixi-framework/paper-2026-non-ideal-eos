const FIGDIR = joinpath(@__DIR__, "figs")
if !isdir(FIGDIR)
    mkpath(FIGDIR)
end

using Plots
using LaTeXStrings
using Trixi

include("trixi_utils.jl")

# compare density profiles for central vs APEC
trixi_include("trixi_transcritical_wave.jl",
            eos_setup = :(:VanDerWaals),
            stabilization = :(:apec),
            initial_refinement_level = 2,
            cfl = 0.5,
            volume_flux = flux_terashima_etal,
            surface_flux = flux_terashima_etal)

x = LinRange(-0.5, 0.5, 1000)
u_exact = @. initial_condition(SVector(x), sol.t[end], equations)
p = plot(x, getindex.(u_exact, 1), label="Exact", linestyle=:solid, linewidth=2, markersize=4, dpi=400)

pd = PlotData1D(sol)
plot!(p, pd["rho"], label="Central", linestyle=:dot, linewidth=2, markersize=4, dpi=400)

trixi_include("trixi_transcritical_wave.jl",
            eos_setup = :(:VanDerWaals),
            stabilization = :(:apec),
            cfl = 0.5,
            initial_refinement_level = 2,
            volume_flux = flux_terashima_etal,
            surface_flux = flux_terashima_etal)
pd = PlotData1D(sol)
plot!(p, pd["rho"], label="APEC", linestyle=:dash, linewidth=2, markersize=4, dpi=400)

plot!(p, legendfontsize=12, guidefontsize=14, tickfontsize=12, titlefontsize=14,
      xlabel=L"x", ylabel=L"Density $\rho$", legend=:bottomright)
savefig(joinpath(FIGDIR, "EPEC_comparison_density_profile.png"))


# calculate errors over time
for initial_refinement_level in (2, 3)

    p1 = plot()
    p2 = plot()
    trixi_include("trixi_transcritical_wave.jl",
                eos_setup = :(:VanDerWaals),
                stabilization = :(:apec),
                initial_refinement_level = :($initial_refinement_level),
                cfl = 0.5,
                volume_flux = flux_central,
                surface_flux = flux_central)
    u, x = sol_and_coordinates(sol)
    u_exact = @. initial_condition(SVector(x), sol.t[end], equations)
    plot!(p1, sol.t, calc_error_1d(sol), label="Central",
        linestyle=:solid, linewidth=2, markersize=4, dpi=400)

    plot!(p2, sol.t, calc_field_error_1d(sol; field=pressure), label="Central",
        linestyle=:solid, linewidth=2, markersize=4, dpi=400, yaxis=:log)

    trixi_include("trixi_transcritical_wave.jl",
                eos_setup = :(:VanDerWaals),
                stabilization = :(:apec),
                initial_refinement_level = :($initial_refinement_level),
                cfl = 0.001,
                volume_flux = flux_coppola,
                surface_flux = flux_coppola)
    u, x = sol_and_coordinates(sol)
    u_exact = @. initial_condition(SVector(x), sol.t[end], equations)
    plot!(p1, sol.t, calc_error_1d(sol), label="EPEC (Coppola et al. 2026)",
        linestyle=:solid, linewidth=2, markersize=4, dpi=400)

    plot!(p2, sol.t, calc_field_error_1d(sol; field=pressure), label="EPEC (Coppola et al. 2026)",
        linestyle=:solid, linewidth=2, markersize=4, dpi=400, yaxis=:log)


    trixi_include("trixi_transcritical_wave.jl",
                eos_setup = :(:VanDerWaals),
                stabilization = :(:apec),
                cfl = 0.5,
                initial_refinement_level = :($initial_refinement_level),
                volume_flux = flux_epec_ranocha,
                surface_flux = flux_epec_ranocha)
    u, x = sol_and_coordinates(sol)
    u_exact = @. initial_condition(SVector(x), sol.t[end], equations)
    plot!(p1, sol.t, calc_error_1d(sol), label="EPEC (correction)",
        linestyle=:dash, linewidth=2, markersize=4, dpi=400)
    plot!(p2, sol.t, calc_field_error_1d(sol; field=pressure), label="EPEC (correction)",
        linestyle=:dash, linewidth=2, markersize=4, dpi=400, yaxis=:log)

    trixi_include("trixi_transcritical_wave.jl",
                eos_setup = :(:VanDerWaals),
                stabilization = :(:apec),
                cfl = 0.5,
                initial_refinement_level = :($initial_refinement_level),
                volume_flux = flux_terashima_etal,
                surface_flux = flux_terashima_etal)
    u, x = sol_and_coordinates(sol)
    u_exact = @. initial_condition(SVector(x), sol.t[end], equations)
    plot!(p1, sol.t, calc_error_1d(sol), label="APEC",
        linestyle=:dot, linewidth=2, markersize=4, dpi=400)
    plot!(p2, sol.t, calc_field_error_1d(sol; field=pressure), label="APEC",
        linestyle=:dot, linewidth=2, markersize=4, dpi=400, yaxis=:log)

    plot!(p1, legendfontsize=12, guidefontsize=14, tickfontsize=12, titlefontsize=14,
        xlabel=L"t", ylabel=L"$L^2$ error", legend=:bottomright)
    plot!(p2, legendfontsize=12, guidefontsize=14, tickfontsize=12, titlefontsize=14,
        xlabel=L"t", ylabel=L"$L^2$ pressure error", legend=:bottomright, ylims=(1e-12, 1e-2))
    savefig(p1, joinpath(FIGDIR, "EPEC_comparison_refinement_level_$(initial_refinement_level).png"))
    savefig(p2, joinpath(FIGDIR, "EPEC_comparison_refinement_level_$(initial_refinement_level)_pressure.png"))

end
