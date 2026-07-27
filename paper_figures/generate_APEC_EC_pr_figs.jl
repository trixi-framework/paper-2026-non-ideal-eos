const FIGDIR = joinpath(@__DIR__, "figs")
if !isdir(FIGDIR)
    mkpath(FIGDIR)
end

using Plots
using Trixi
using LaTeXStrings

include("trixi_utils.jl")

for initial_condition in (initial_condition_transcritical_wave, initial_condition_transcritical_wave_sharp)
  for (polydeg, initial_refinement_level) in ((3, 4), (7, 3),)

      p1 = plot()
      trixi_include("trixi_transcritical_wave.jl", 
                  eos_setup = :(:PengRobinson),
                  initial_condition = initial_condition,
                  stabilization = :(:apec),
                  cfl = 0.5,
                  polydeg = polydeg,
                  initial_refinement_level = :($initial_refinement_level),
                  volume_flux = flux_terashima_etal, 
                  surface_flux = FluxPlusDissipation(flux_terashima_etal, DissipationLocalLaxFriedrichs()), 
                  tspan = (0.0, 0.1))
      u, x = sol_and_coordinates(sol)
      u_exact = @. initial_condition(SVector(x), sol.t[end], equations)
      plot!(p1, sol.t, calc_error_1d(sol), label="APEC", 
          linestyle=:dot, linewidth=2, markersize=4, dpi=400)

      trixi_include("trixi_transcritical_wave.jl", 
                  eos_setup = :(:PengRobinson),
                  initial_condition = initial_condition,
                  stabilization = :(:apec_ecav),
                  cfl = 0.5,
                  polydeg = polydeg,
                  initial_refinement_level = :($initial_refinement_level),
                  volume_flux = flux_terashima_etal, 
                  surface_flux = FluxPlusDissipation(flux_terashima_etal, DissipationLocalLaxFriedrichs()),
                  tspan = (0.0, 0.1))
      u, x = sol_and_coordinates(sol)
      u_exact = @. initial_condition(SVector(x), sol.t[end], equations)
      plot!(p1, sol.t, calc_error_1d(sol), label="APEC + EC", 
          linestyle=:dashdot, linewidth=2, markersize=4, dpi=400)
      plot!(p1, max_alpha_cb.affect!.times, max.(eps(), max_alpha_cb.affect!.max_alpha), 
            label="Blending coefficient", yaxis=:log)


      if initial_condition == initial_condition_transcritical_wave
        plot!(p1, legendfontsize=12, guidefontsize=14, tickfontsize=12, titlefontsize=14, 
              xlabel=L"t", ylabel=L"$L^2$ error", legend=:topleft, yaxis=:log, ylims=(1e-9, 1e-1))
        savefig(p1, joinpath(FIGDIR, "blending_coeff_smoothwave_polydeg_$(polydeg)_refinement_level_$(initial_refinement_level).png"))
      else
        plot!(p1, legendfontsize=12, guidefontsize=14, tickfontsize=12, titlefontsize=14, 
              xlabel=L"t", ylabel=L"$L^2$ error", legend=:topleft, yaxis=:log, ylims=(1e-4, 1e1))
        savefig(p1, joinpath(FIGDIR, "blending_coeff_sharpwave_polydeg_$(polydeg)_refinement_level_$(initial_refinement_level).png"))
      end

  end

end
