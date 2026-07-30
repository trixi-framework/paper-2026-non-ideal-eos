const FIGDIR = joinpath(@__DIR__, "figs")
if !isdir(FIGDIR)
    mkpath(FIGDIR)
end

using Plots
using LaTeXStrings
using LinearAlgebra: eigvals
using Trixi
using Trixi: ForwardDiff

include("trixi_utils.jl")

# `Trixi.jacobian_ad_forward` differentiates with an untagged `ForwardDiff.JacobianConfig`.
# The thermodynamic routines (e.g. `heat_capacity_constant_volume`) differentiate the equation
# of state with `ForwardDiff.derivative` themselves, and the resulting nested dual numbers
# cannot be ordered against untagged ones. Hence, we redo the same Jacobian with a tag.
struct SemiJacobianTag end

function jacobian_ad(semi, t0, u0_ode)
    du_ode = similar(u0_ode)
    config = ForwardDiff.JacobianConfig(SemiJacobianTag(), du_ode, u0_ode)
    semi_dual = Trixi.remake(semi, uEltype = eltype(config))
    # the tag of `config` does not belong to the closure below, so skip the tag check
    return ForwardDiff.jacobian(du_ode, u0_ode, config, Val{false}()) do du_ode, u_ode
        return Trixi.rhs!(du_ode, u_ode, semi_dual, t0)
    end
end

# eigenvalues of the Jacobian of the semidiscretization along the solution trajectory.
# The snapshots are subsampled with a fixed stride (instead of a fixed number of frames)
# such that the snapshot times of different runs coincide, also if a run stopped early.
function eigenvalue_spectra(sol, semi; nframes=100)
    indices = 1:max(1, cld(length(sol.u), nframes)):length(sol.u)
    lambda = map(indices) do i
        return eigvals(jacobian_ad(semi, sol.t[i], sol.u[i]))
    end
    return (; t = sol.t[indices], lambda)
end

# animate the spectra of several runs next to each other, where `spectra` is a vector of
# `label => eigenvalue_spectra(sol, semi)` pairs. Each panel keeps its own axis limits fixed
# over time so that the frames are comparable. The limits differ between the panels since the
# spectra of the different fluxes differ by orders of magnitude.
function animate_spectra(filename, spectra; fps=10, layout=(2, 2), size=(1000, 800))
    function lims(values)
        lo, hi = extrema(values)
        pad = 0.05 * (hi - lo) + eps(max(abs(lo), abs(hi)))
        return (lo - pad, hi + pad)
    end
    xlims = [lims(real.(reduce(vcat, s.lambda))) for (_, s) in spectra]
    ylims = [lims(imag.(reduce(vcat, s.lambda))) for (_, s) in spectra]

    # the runs may have different lengths if one of them blew up before `tspan[end]`
    nframes = maximum(length(s.lambda) for (_, s) in spectra)
    frame_times = reduce((t1, t2) -> length(t1) >= length(t2) ? t1 : t2, s.t for (_, s) in spectra)

    anim = @animate for frame in 1:nframes
        panels = map(enumerate(spectra)) do (j, (label, s))
            i = min(frame, length(s.lambda))
            lambda = s.lambda[i]
            title = i == frame ? label : label * " (stopped at t = $(round(s.t[i], sigdigits=3)))"
            p = scatter(real.(lambda), imag.(lambda), markersize=3, markerstrokewidth=0,
                        label=L"\max \mathrm{Re}\, \lambda = %$(round(maximum(real, lambda), sigdigits=3))",
                        title=title, xlims=xlims[j], ylims=ylims[j], legend=:topleft,
                        xlabel=L"\mathrm{Re}\, \lambda", ylabel=L"\mathrm{Im}\, \lambda",
                        # rotated ticks since the real parts can span several orders of magnitude
                        xrotation=20,
                        legendfontsize=10, guidefontsize=12, tickfontsize=10, titlefontsize=12)
            vline!(p, [0], color=:black, linestyle=:dash, label=nothing)
            return p
        end
        plot(panels..., layout=layout, size=size,
             plot_title=L"t = %$(round(frame_times[frame], sigdigits=3))")
    end
    return gif(anim, filename, fps=fps)
end

# plot a scalar `measure` of the spectrum (e.g. `-min Re λ` or `max |λ|`) over time for all
# runs in `spectra`, where `measure` maps a vector of eigenvalues to a positive number
function plot_spectrum_measure(filename, spectra, measure, ylabel;
                               linestyles=(:solid, :solid, :dash, :dot))
    p = plot()
    values = [[measure(lambda) for lambda in s.lambda] for (_, s) in spectra]
    for (j, (label, s)) in enumerate(spectra)
        plot!(p, s.t, values[j], label=label,
              linestyle=linestyles[mod1(j, length(linestyles))],
              linewidth=2, markersize=4, dpi=400, yaxis=:log)
    end

    # the central flux has an (up to round-off) purely imaginary spectrum for the symmetric
    # initial condition, so `-min Re λ` would stretch the logarithmic axis over ~20 decades.
    # Hence, everything below eight decades of the largest value is cut off.
    upper = maximum(maximum, values)
    lower = max(minimum(minimum, values), 1e-8 * upper)

    # the strongly oscillating curves leave no free space for the legend inside the axes
    plot!(p, legendfontsize=12, guidefontsize=14, tickfontsize=12, titlefontsize=14,
          xlabel=L"t", ylabel=ylabel, legend=:outertop, legend_columns=2,
          ylims=(0.5 * lower, 2 * upper), size=(800, 500))
    savefig(p, filename)
    return p
end

# the most negative real part of the spectrum is the fastest decaying mode, i.e. it sets the
# explicit time step and thus measures the stiffness of the semidiscretization. It is negative
# for all fluxes considered here, so `-min Re λ` is positive.
plot_stiffness(filename, spectra) = plot_spectrum_measure(filename, spectra,
                                                          lambda -> -minimum(real, lambda),
                                                          L"-\min \mathrm{Re}\, \lambda")

# the spectral radius bounds the time step of the explicit Runge-Kutta method
plot_spectral_radius(filename, spectra) = plot_spectrum_measure(filename, spectra,
                                                                lambda -> maximum(abs, lambda),
                                                                L"\max |\lambda|")

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
    spectra_central = eigenvalue_spectra(sol, semi)

    trixi_include("trixi_transcritical_wave.jl",
                eos_setup = :(:VanDerWaals),
                stabilization = :(:apec),
                initial_refinement_level = :($initial_refinement_level),
                cfl = 0.01,
                volume_flux = flux_coppola,
                surface_flux = flux_coppola)
    u, x = sol_and_coordinates(sol)
    u_exact = @. initial_condition(SVector(x), sol.t[end], equations)
    plot!(p1, sol.t, calc_error_1d(sol), label="EPEC (Coppola et al. 2026)",
        linestyle=:solid, linewidth=2, markersize=4, dpi=400)

    plot!(p2, sol.t, calc_field_error_1d(sol; field=pressure), label="EPEC (Coppola et al. 2026)",
        linestyle=:solid, linewidth=2, markersize=4, dpi=400, yaxis=:log)
    spectra_coppola = eigenvalue_spectra(sol, semi)


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
    spectra_epec_ranocha = eigenvalue_spectra(sol, semi)

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
    spectra_apec = eigenvalue_spectra(sol, semi)

    spectra = ["Central" => spectra_central,
               "EPEC (Coppola et al. 2026)" => spectra_coppola,
               "EPEC (correction)" => spectra_epec_ranocha,
               "APEC" => spectra_apec]
    animate_spectra(joinpath(FIGDIR, "EPEC_spectra_refinement_level_$(initial_refinement_level).gif"),
                    spectra)
    plot_stiffness(joinpath(FIGDIR, "EPEC_stiffness_refinement_level_$(initial_refinement_level).png"),
                   spectra)
    plot_spectral_radius(joinpath(FIGDIR, "EPEC_spectral_radius_refinement_level_$(initial_refinement_level).png"),
                         spectra)

    plot!(p1, legendfontsize=12, guidefontsize=14, tickfontsize=12, titlefontsize=14,
        xlabel=L"t", ylabel=L"$L^2$ error", legend=:bottomright)
    plot!(p2, legendfontsize=12, guidefontsize=14, tickfontsize=12, titlefontsize=14,
        xlabel=L"t", ylabel=L"$L^2$ pressure error", legend=:bottomright, ylims=(1e-15, 1e-2))
    savefig(p1, joinpath(FIGDIR, "EPEC_comparison_refinement_level_$(initial_refinement_level).png"))
    savefig(p2, joinpath(FIGDIR, "EPEC_comparison_refinement_level_$(initial_refinement_level)_pressure.png"))

end
