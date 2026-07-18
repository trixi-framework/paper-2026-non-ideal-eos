using Plots
using Trixi

include("trixi_utils.jl")

initial_condition = initial_condition_transcritical_wave

for (polydeg, initial_refinement_level) in ((1, 5), (3, 4), (7, 3))

    p1 = plot()
    p2 = plot()
    trixi_include("trixi_transcritical_wave.jl", 
                eos_setup = :(:PengRobinson),
                initial_condition = initial_condition,
                stabilization = :(:apec),
                cfl = 0.5,
                polydeg = polydeg,
                initial_refinement_level = :($initial_refinement_level),
                volume_flux = flux_terashima_etal, 
                surface_flux = flux_terashima_etal, 
                tspan = (0.0, 0.1))
    u, x = sol_and_coordinates(sol)
    u_exact = @. initial_condition(SVector(x), sol.t[end], equations)
    plot!(p1, sol.t, calc_error_1d(sol), label="APEC", 
        linestyle=:dot, linewidth=2, markersize=4, dpi=400)
    plot!(p2, sol.t, calc_field_error_1d(sol; field=pressure), label="APEC", 
        linestyle=:dot, linewidth=2, markersize=4, dpi=400)

    trixi_include("trixi_transcritical_wave.jl", 
                eos_setup = :(:PengRobinson),
                stabilization = :(:apec),
                cfl = 0.5,
                polydeg = polydeg,
                initial_refinement_level = :($initial_refinement_level),
                volume_flux = flux_terashima_etal, 
                surface_flux = FluxPlusDissipation(flux_terashima_etal, DissipationLocalLaxFriedrichs()),
                tspan = (0.0, 0.1))
    u, x = sol_and_coordinates(sol)
    u_exact = @. initial_condition(SVector(x), sol.t[end], equations)
    plot!(p1, sol.t, calc_error_1d(sol), label="APEC + LxF penalty", 
        linestyle=:dashdot, linewidth=2, markersize=4, dpi=400)
    plot!(p2, sol.t, calc_field_error_1d(sol; field=pressure), label="APEC + LxF penalty", 
        linestyle=:dashdot, linewidth=2, markersize=4, dpi=400)

    trixi_include("trixi_transcritical_wave.jl", 
                eos_setup = :(:PengRobinson),
                stabilization = :(:apec),
                cfl = 0.5,
                polydeg = polydeg,
                initial_refinement_level = :($initial_refinement_level),
                volume_flux = flux_terashima_etal, 
                surface_flux = flux_lax_friedrichs,
                tspan = (0.0, 0.1))
    u, x = sol_and_coordinates(sol)
    u_exact = @. initial_condition(SVector(x), sol.t[end], equations)
    plot!(p1, sol.t, calc_error_1d(sol), label="APEC + LxF interface", 
        linestyle=:dash, linewidth=2, markersize=4, dpi=400)
    plot!(p2, sol.t, calc_field_error_1d(sol; field=pressure), label="APEC + LxF interface", 
        linestyle=:dash, linewidth=2, markersize=4, dpi=400)

    using LaTeXStrings
    plot!(p1, legendfontsize=12, guidefontsize=14, tickfontsize=12, titlefontsize=14, 
        xlabel=L"t", ylabel=L"$L^2$ error", legend=:topleft)
    plot!(p2, legendfontsize=12, guidefontsize=14, tickfontsize=12, titlefontsize=14, 
        xlabel=L"t", ylabel=L"$L^2$ pressure error", yaxis=:log, ylims = (1e-8, 1e-0), 
        legend=:bottomright)
    savefig(p1, "figs/APEC_LxF_comparison_polydeg_$(polydeg)_refinement_level_$(initial_refinement_level).png")
    savefig(p2, "figs/APEC_LxF_comparison_pressure_polydeg_$(polydeg)_refinement_level_$(initial_refinement_level).png")
end


