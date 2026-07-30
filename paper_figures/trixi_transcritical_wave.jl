using OrdinaryDiffEqLowStorageRK
using SciMLBase: DiscreteCallback
using Trixi
using Trixi: ForwardDiff

include("trixi_utils.jl")

###############################################################################
# semidiscretization of the compressible Euler equations

# the smooth Peng-Robinson N2 transcritical wave taken from "An entropy-stable hybrid
# scheme for simulations of transcritical real-fluid flows" by Ma, Ihme (2017). In this
# context, the wave is "transcritical" because the solution involves both subcritical
# and supercritical density and temperature values.
#
# <https://doi.org/10.1016/j.jcp.2017.03.022>
function initial_condition_transcritical_wave(x, t,
                                              equations::NonIdealCompressibleEulerEquations1D{<:PengRobinson})
    RealT = eltype(x)
    eos = equations.equation_of_state

    k_wave = 2
    rho_min, rho_max = 56.9, 793.1
    # rho_min, rho_max = 700.0, 793.1
    v1 = 100
    rho = 0.5f0 * (rho_min + rho_max) +
          0.5f0 * (rho_max - rho_min) * sin(k_wave * pi * (x[1] - v1 * t))
    p = 5e6 # transcritical
    # p = 5e7 # supercritical

    V = inv(rho)

    # invert for temperature given p, V
    T = temperature_pV(p, V, eos)

    return thermo2cons(SVector(V, v1, T), equations)
end

function initial_condition_transcritical_wave_sharp(x, t,
                                              equations::NonIdealCompressibleEulerEquations1D)
    RealT = eltype(x)
    eos = equations.equation_of_state

    rho_min, rho_max = 56.9, 793.1
    v1 = 100
    rho = 0.5 * (rho_min + rho_max) +
        0.5 * (rho_max - rho_min) * tanh(20 * (abs(sin(pi * (x[1] - v1 * t))) - 0.75))
    p = 5e6

    V = inv(rho)

    # invert for temperature given p, V
    T = temperature_pV(p, V, eos)

    return thermo2cons(SVector(V, v1, T), equations)
end

function initial_condition_density_wave_coppola(x, t, equations::NonIdealCompressibleEulerEquations1D{<:VanDerWaals})
    eos = equations.equation_of_state

    A, B = 0.07, 0.12
    rho_0 = inv(3 * eos.b) # critical density
    v1 = 1.0
    rho = rho_0 * (A + B * exp(sin(2 * pi * (x[1] - v1 * t))))
    p = 100 # supercritical

    V = inv(rho)

    # invert for temperature given p, V
    T = temperature_pV(p, V, eos)

    return thermo2cons(SVector(V, v1, T), equations)
end

###########################

eos_setup = :PengRobinson
if eos_setup == :PengRobinson
    eos = PengRobinson()
    initial_condition = initial_condition_transcritical_wave
    # initial_condition = initial_condition_transcritical_wave_sharp
    tspan = (0.0, 0.1)
else
    eos = VanDerWaals(; a = 5.94768233178e-3, b = 1.72768204288e-3, gamma = 1.4, R = 1) # from Coppola
    initial_condition = initial_condition_density_wave_coppola
    tspan = (0.0, 2.0)
end

equations = NonIdealCompressibleEulerEquations1D(eos)

# volume_flux = flux_central
# volume_flux = flux_coppola
# volume_flux = flux_epec_ranocha
volume_flux = flux_terashima_etal

# surface_flux = flux_terashima_etal
# surface_flux = flux_central
# surface_flux = flux_lax_friedrichs
# surface_flux = FluxPlusDissipation(volume_flux, dissipation_barth)
# surface_flux = FluxPlusDissipation(volume_flux, DissipationLocalLaxFriedrichs())
surface_flux = volume_flux

polydeg = 3
basis = LobattoLegendreBasis(polydeg)
coordinates_min, coordinates_max = -0.5, 0.5
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 4,
                periodicity = true)

stabilization = :apec
stabilization = :apec_ecav
# stabilization = :shock_capturing

###########################

# pure APEC
volume_integral_apec = VolumeIntegralFluxDifferencing(volume_flux)

# APEC + ECAV
indicator = IndicatorEntropyCorrection(equations, basis)
volume_integral_default = volume_integral_apec
volume_integral_entropy_stable = VolumeIntegralPureLGLFiniteVolume(surface_flux)
volume_integral_apec_ecav = VolumeIntegralAdaptive(indicator,
                                                   volume_integral_default,
                                                   volume_integral_entropy_stable)

if stabilization == :apec
    volume_integral = volume_integral_apec
elseif stabilization == :apec_ecav
    volume_integral = volume_integral_apec_ecav
elseif stabilization == :shock_capturing
    indicator_sc = IndicatorHennemannGassner(equations, basis,
                                             alpha_max = 0.05,
                                             alpha_min = 0.00,
                                             alpha_smooth = true,
                                             variable = Trixi.density_pressure)
    volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                     volume_flux_dg = volume_flux,
                                                     volume_flux_fv = surface_flux)
end

dg = DGSEM(basis, surface_flux, volume_integral)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, dg;
                                    boundary_conditions = boundary_condition_periodic)

###############################################################################
# max-alpha callback (ECAV / shock capturing only)

struct MaxAlphaCallback
    indicator
    semi
    times::Vector{Float64}
    max_alpha::Vector{Float64}
    MaxAlphaCallback(indicator, semi) = new(indicator, semi, Float64[], Float64[])
end

function (cb::MaxAlphaCallback)(integrator)
    alpha = cb.indicator.cache.alpha
    isempty(alpha) && return nothing
    solution_norm = calc_solution_norm_1d(integrator.u, cb.semi)
    push!(cb.times, integrator.t)
    # push!(cb.max_alpha, maximum(alpha) / solution_norm)
    push!(cb.max_alpha, maximum(alpha))
    SciMLBase.u_modified!(integrator, false)
    return nothing
end

function max_alpha_callback(indicator, semi)
    condition = (u, t, integrator) -> true
    return DiscreteCallback(condition, MaxAlphaCallback(indicator, semi),
                            save_positions = (false, false))
end

###############################################################################
# ODE solvers, callbacks etc.

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()
alive_callback = AliveCallback(alive_interval = 5000)
analysis_callback = AnalysisCallback(semi, interval = 20000)

###############################################################################
# run the simulation

stepsize_callback = StepsizeCallback(cfl = 0.5)
solver = CarpenterKennedy2N54(williamson_condition = false)
if stabilization == :apec_ecav
    max_alpha_cb = max_alpha_callback(indicator, semi)
    callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback,
                            stepsize_callback, max_alpha_cb)
elseif stabilization == :shock_capturing
    max_alpha_cb = max_alpha_callback(indicator_sc, semi)
    callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback,
                            stepsize_callback, max_alpha_cb)
else
    max_alpha_cb = nothing
    callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback,
                            stepsize_callback)
end
kwargs = (; maxiters=1_000_000)

sol = solve(ode, solver; dt = stepsize_callback(ode), # solve needs some value here but it will be overwritten by the stepsize_callback
            kwargs...,
            ode_default_options()..., saveat=LinRange(tspan..., 500),
            callback = callbacks);

if sol.retcode !== ReturnCode.Success
    println("sol.retcode is not successful = $(sol.retcode)")
end

# using Plots
# plot(sol.t, calc_error_1d(sol))

# if eos isa PengRobinson
#     p1 = plot(PlotData1D(sol)["rho"], label="p=$polydeg", linewidth=2)
#     plot!(p1, vec(x), vec(getindex.(u_exact, 1)), label="Exact", linewidth=2)
#     p2 = plot(PlotData1D(sol)["p"], ylims=(1e4, 1e7), label="p=$polydeg", linewidth=2)
#     plot!(p2, SVector(extrema(x)), 5e6 * [1;1], label="Exact", linewidth=2)
#     plot(p1, p2, leg=true)
#     # plot(sol)
# else
#     p1 = plot(PlotData1D(sol)["rho"], label="p=$polydeg", linewidth=2)
#     plot!(p1, vec(x), vec(getindex.(u_exact, 1)), label="Exact", linewidth=2)
#     p2 = plot(PlotData1D(sol)["p"], label="p=$polydeg", linewidth=2)
#     plot!(p1, p2, leg=true)
# end

# @gif for (i, u) in enumerate(sol.u)
#     u_exact = @. initial_condition(SVector(x), sol.t[i], equations)
#     p1 = plot(PlotData1D(u, semi)["rho"], label="p=$polydeg", linewidth=2)
#     plot!(p1, vec(x), vec(getindex.(u_exact, 1)), label="Exact", linewidth=2, legend=:topright)
#     p2 = plot(PlotData1D(u, semi)["p"], ylims=(1e6, 1e7), label="p=$polydeg", linewidth=2)
#     plot!(p2, SVector(extrema(x)), 5e6 * [1;1], label="Exact", linewidth=2)
#     plot(p1, p2, leg=true, legend=:topright)
# end
