using OrdinaryDiffEqSSPRK, OrdinaryDiffEqLowStorageRK
using Trixi
using Trixi: ForwardDiff

###############################################################################
# semidiscretization of the compressible Euler equations

eos = PengRobinson()
equations = NonIdealCompressibleEulerEquations2D(eos)

function initial_condition_jet(x, t, equations::NonIdealCompressibleEulerEquations2D)
    eos = equations.equation_of_state

    v1, v2 = 0.0, 0.0

    rho = 45 # kg/m3
    T = 290.2 # K, about 4 MPa 

    # # from Ma et al
    # rho = 45.5
    # T = 300

    V = inv(rho)
    e_internal = energy_internal_specific(V, T, eos)

    return SVector(rho, rho * v1, rho * v2, rho * (e_internal + 0.5 * (v1^2 + v2^2)))
end

function boundary_condition_jet(x, t, equations::NonIdealCompressibleEulerEquations2D)
    eos = equations.equation_of_state
    h = 0.0022 # 2.2 mm
    
    # Smoothing this discontinuous initial profile based on grid spacing
    # add a small v2 perturbation to break symmetry
    v2_perturbation = 1e-4
    q_bulk = SVector(500, 100, v2_perturbation, 124.6) # rho, v1, v2, T
    q_in = SVector(45, 0.0, 0.0, 290.2)

    function prim2cons_jet(q, eos)
        rho, v1, v2, T = q
        V = inv(rho)
        e_internal = energy_internal_specific(V, T, eos)
        return SVector(rho, rho * v1, rho * v2, rho * (e_internal + 0.5 * (v1^2 + v2^2)))
    end

    smooth_step(y) = 0.5 * (tanh(500 * (y + 0.5 * h) / h) - tanh(500 * (y - 0.5 * h) / h))
    return (1 - smooth_step(x[2])) * prim2cons_jet(q_in, eos) + smooth_step(x[2]) * prim2cons_jet(q_bulk, eos)
end

# Newton's relative residual seems to stall at 5e-8 for this problem
# and at ~6 iterations. 
@inline Trixi.eos_newton_tol(eos::PengRobinson) = 6e-8
@inline Trixi.eos_newton_maxiter(eos::PengRobinson) = 10
@inline Trixi.eos_initial_temperature(V, e_internal, eos::PengRobinson) = eos.Tc

initial_condition = initial_condition_jet

volume_flux = flux_terashima_etal
surface_flux = FluxPlusDissipation(volume_flux, DissipationLocalLaxFriedrichs())
basis = LobattoLegendreBasis(3)

indicator_ec = IndicatorEntropyCorrection(equations, basis; 
                                          scaling = 1.5, 
                                          alpha_smooth=true)

# @inline function inv_density_temperature(u, equations::NonIdealCompressibleEulerEquations2D{<:PengRobinson})
#     eos = equations.equation_of_state
#     V, _, _, T = cons2thermo(u, equations)
#     # flags whenver V ≈ b and T ≈ 0
#     return inv((V - eos.b) * T) 
# end
# 
# indicator_sc = IndicatorHennemannGassner(equations, basis,
#                                          alpha_max = 0.01,
#                                          alpha_min = 0.0,
#                                          alpha_smooth = false,
#                                          variable = inv_density_temperature)
# indicator = IndicatorEntropyCorrectionShockCapturingCombined(; indicator_entropy_correction=indicator_ec,
#                                                                indicator_shock_capturing=indicator_sc)
indicator = indicator_ec
volume_integral_apec = VolumeIntegralFluxDifferencing(volume_flux)
volume_integral_entropy_stable = VolumeIntegralPureLGLFiniteVolume(surface_flux)
volume_integral = VolumeIntegralAdaptive(indicator,
                                         volume_integral_apec,
                                         volume_integral_entropy_stable)

solver = DGSEM(basis, surface_flux, volume_integral)

h = 0.0022 # 2.2 mm
coordinates_min = (0.0, -8 * h)
coordinates_max = (32 * h, 8 * h)
cells_per_dimension = (320, 160)
mesh = StructuredMesh(cells_per_dimension,
                      coordinates_min, coordinates_max,
                      periodicity = (false, true))

boundary_condition_inflow = BoundaryConditionDirichlet(boundary_condition_jet)
boundary_condition_outflow = boundary_condition_slip_wall

boundary_conditions = (x_neg = boundary_condition_inflow,
                       x_pos = boundary_condition_outflow,
                       y_neg = boundary_condition_periodic,
                       y_pos = boundary_condition_periodic)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver;
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 1e-3)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 5000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

function gamma(u, equations::NonIdealCompressibleEulerEquations2D{<:PengRobinson})
    eos = equations.equation_of_state
    V, _, _, T = cons2thermo(u, equations)

    (; cv0, R) = eos
    dpdT_V, dpdV_T = Trixi.calc_pressure_derivatives(V, T, eos)
    K1   = Trixi.calc_K1(V, eos)
    d2aT = Trixi.peng_robinson_d2a(T, eos)

    cv = cv0 - K1 * T * d2aT
    cp = (cv0 + R) - R - K1 * T * d2aT - T * dpdT_V^2 / dpdV_T

    return cp / cv
end

function Trixi.get_node_variable(::Val{:gamma}, u, mesh,
                                 equations::NonIdealCompressibleEulerEquations2D{<:PengRobinson},
                                 dg, cache)
    n_nodes = nnodes(dg)
    n_elements = nelements(dg, cache)
    gamma_array = zeros(eltype(cache.elements), n_nodes, n_nodes, n_elements)

    Trixi.@threaded for element in eachelement(dg, cache)
        for j in eachnode(dg), i in eachnode(dg)
            u_node = get_node_vars(u, equations, dg, i, j, element)
            gamma_array[i, j, element] = gamma(u_node, equations)
        end
    end

    return gamma_array
end

output_directory = "out"
save_solution = SaveSolutionCallback(interval = 10000,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim,
                                     extra_node_variables = (:gamma,),
                                     output_directory = output_directory)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

stepsize_callback = StepsizeCallback(cfl = 0.8)

###############################################################################
# run the simulation

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback, 
                        stepsize_callback, 
                        save_solution)

ode_solver = CarpenterKennedy2N54(williamson_condition = false)
sol = solve(ode, ode_solver;
            adaptive=false, dt = stepsize_callback(ode), maxiters = 5_000_000,
            ode_default_options()..., callback = callbacks);

using Plots
plot(ScalarPlotData2D(Trixi.density, sol.u[end], semi), ylims=(-8 * h, 8 * h), dpi = 400)
savefig("polydeg$(Trixi.polydeg(basis))_transcritical_jet_rho.png")

# x = semi.cache.elements.node_coordinates[1, ..]
# y = semi.cache.elements.node_coordinates[2, ..]
# x = reshape(x, :, Trixi.nelements(mesh, solver, semi.cache))
# y = reshape(y, :, Trixi.nelements(mesh, solver, semi.cache))
# xc = vec(sum(x, dims = 1) / size(x, 2))
# yc = vec(sum(y, dims = 1) / size(y, 2))
# alpha = copy(indicator_ec.cache.alpha)
# scatter(xc, yc, zcolor=alpha, msw=0, ms = 4)

using Trixi2Vtk
trixi2vtk(joinpath(output_directory, "solution_*.h5");
          output_directory = joinpath(output_directory, "vtk"),
          verbose = true)
