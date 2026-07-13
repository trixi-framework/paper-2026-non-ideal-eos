using OrdinaryDiffEqLowStorageRK
using Trixi
using Trixi: ForwardDiff

###############################################################################
# semidiscretization of the compressible Euler equations

eos = PengRobinson()
equations = NonIdealCompressibleEulerEquations1D(eos)

# the Peng-Robinson N2 transcritical shock taken from "An entropy-stable hybrid scheme 
# for simulations of transcritical real-fluid flows" by Ma, Ihme (2017).
# <https://doi.org/10.1016/j.jcp.2017.03.022>
function initial_condition_transcritical_shock(x, t,
                                               equations::NonIdealCompressibleEulerEquations1D{<:PengRobinson})
    RealT = eltype(x)
    eos = equations.equation_of_state

    u_ll = SVector(800, 0, convert(RealT, 60.0e6))
    u_rr = SVector(80, 0, convert(RealT, 6.0e6))
    if x[1] < 0
        rho, v1, p = u_ll
    elseif x[1] ≈ 0
        rho, v1, p = 0.5 * (u_ll + u_rr)
    else
        rho, v1, p = u_rr
    end 

    V = inv(rho)

    # invert for temperature given p, V
    T = temperature_given_Vp(V, p, eos)

    return thermo2cons(SVector(V, v1, T), equations)
end
initial_condition = initial_condition_transcritical_shock

volume_flux = flux_central
volume_flux = flux_central_terashima_etal
volume_flux = flux_terashima_etal

surface_flux = FluxPlusDissipation(volume_flux, DissipationLocalLaxFriedrichs())

polydeg = 3
initial_refinement_level = 7
polydeg, initial_refinement_level = 7, 6

basis = LobattoLegendreBasis(polydeg)
indicator_ec = IndicatorEntropyCorrection(equations, basis; 
                                          scaling = true, 
                                          alpha_smooth = false)
volume_integral_apec = VolumeIntegralFluxDifferencing(volume_flux)
volume_integral_entropy_stable = VolumeIntegralPureLGLFiniteVolume(surface_flux)
volume_integral = VolumeIntegralAdaptive(indicator_ec, 
                                         volume_integral_apec, 
                                         volume_integral_entropy_stable)

dg = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = -0.5
coordinates_max = 0.5
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = initial_refinement_level,
                n_cells_max = 30_000,
                periodicity = false)

boundary_conditions = (x_neg = BoundaryConditionDirichlet(initial_condition),
                       x_pos = BoundaryConditionDirichlet(initial_condition))

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, dg,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 5.0e-4)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 2000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 100,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

stepsize_callback = StepsizeCallback(cfl = 0.3)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        stepsize_callback)

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = stepsize_callback(ode), # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks);
if !(sol.retcode == ReturnCode.Success)
    @warn "Failure with return code $(sol.retcode)"
end
