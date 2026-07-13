using OrdinaryDiffEqLowStorageRK
using OrdinaryDiffEqSSPRK
using Trixi
using Trixi: ForwardDiff

###############################################################################
# semidiscretization of the compressible Euler equations

eos = PengRobinson()
equations = NonIdealCompressibleEulerEquations2D(eos)

# A transcritical mixing layer adapted from "Kinetic-energy- and pressure-equilibrium-preserving 
# schemes for real-gas turbulence in the transcritical regime" by Bernades, Jofre, Capuano (2023). 
# <https://doi.org/10.1016/j.jcp.2023.112477>
function initial_condition_transcritical_mixing(x, t,
                                                equations::NonIdealCompressibleEulerEquations2D)
    eos = equations.equation_of_state

    RealT = eltype(x)

    x, y = x

    pc = 3.4e6 # critical pressure for N2
    p = 2 * pc

    # from Bernades et al
    epsilon, delta, A = 1, convert(RealT, 1 / 20), convert(RealT, 3 / 8)
    u0 = 25 # m/s
    T = eos.Tc * (3 * A - A * tanh(y / delta)) # Tc is 126.2 for N2

    tol = Trixi.eos_newton_tol(eos)

    # invert for V given p, T. Initialize V so that the denominator 
    # (V - b) in Peng-Robinson is positive. 
    V = eos.b + tol
    dp = pressure(V, T, eos) - p
    iter = 1
    while abs(dp) > tol * abs(p) && iter < 100
        dp = pressure(V, T, eos) - p
        dpdV_T = ForwardDiff.derivative(V -> pressure(V, T, eos), V)
        V = max(eos.b + tol, V - dp / dpdV_T)
        iter += 1
    end
    if iter == 100
        @warn "Solver for temperature(V, p) did not converge"
    end

    k = 6
    dv = epsilon * sinpi(k * x) * (tanh(100 * (y + 0.1)) - tanh(100 * (y - 0.1))) / 2
    v1 = u0 * (1 + convert(RealT, 0.2) * tanh(y / delta)) + dv
    v2 = dv

    return thermo2cons(SVector(V, v1, v2, T), equations)
end

initial_condition = initial_condition_transcritical_mixing

Trixi.eos_newton_tol(eos) = 1e-8

volume_flux = flux_terashima_etal
surface_flux = FluxPlusDissipation(volume_flux, DissipationLocalLaxFriedrichs())

polydeg = 7
basis = LobattoLegendreBasis(polydeg)
volume_integral_apec = VolumeIntegralFluxDifferencing(volume_flux)
volume_integral_entropy_stable = VolumeIntegralPureLGLFiniteVolume(surface_flux)
indicator = IndicatorEntropyCorrection(equations, basis)
volume_integral = VolumeIntegralAdaptive(indicator,
                                         volume_integral_apec,
                                         volume_integral_entropy_stable)

# volume_integral = volume_integral_default # APEC only
solver = DGSEM(basis, surface_flux, volume_integral)

if Trixi.polydeg(basis) == 3
    cells_per_dimension = (64, 32)
else # if polydeg == 7
    cells_per_dimension = (32, 16)
end
coordinates_min = (-0.5, -0.25)
coordinates_max = (0.5, 0.25)
mesh = StructuredMesh(cells_per_dimension,
                      coordinates_min, coordinates_max,
                      periodicity = (true, false))

boundary_conditions = (x_neg = boundary_condition_periodic,
                       x_pos = boundary_condition_periodic,
                       y_neg = boundary_condition_slip_wall,
                       y_pos = boundary_condition_slip_wall)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0, 0.033 * 2)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

if Trixi.polydeg(basis) == 3
    stepsize_callback = StepsizeCallback(cfl = 0.9)
else # if polydeg == 7
    stepsize_callback = StepsizeCallback(cfl = 0.6)
end

callbacks = CallbackSet(summary_callback,
                        analysis_callback, 
                        alive_callback,
                        stepsize_callback,
                        )

###############################################################################
# run the simulation

solver = CarpenterKennedy2N54(williamson_condition = false)
sol = solve(ode, solver;
            dt = stepsize_callback(ode), # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks);

save_plots = true
if save_plots
    # Evaluate pressure at the PR critical state (Zc ≈ 0.307)
    using Trixi: ForwardDiff
    function critical_compressibility(eos::PengRobinson)
        (; R, b, Tc) = eos
        dpdV(V) = ForwardDiff.derivative(V -> Trixi.pressure(V, Tc, eos), V)

        # Newton solve for Vc with (∂p/∂V)|_{T=Tc} = 0
        Vc = 3.95 * b   # good initial guess for standard PR
        for _ in 1:20
            f  = dpdV(Vc)
            df = ForwardDiff.derivative(dpdV, Vc)
            Vc -= f / df
        end

        Pc = Trixi.pressure(Vc, Tc, eos)
        return Pc * Vc / (R * Tc), Vc, Pc
    end

    # calculate relative velocity and pressure
    function cons2relative(u, equations::NonIdealCompressibleEulerEquations2D)
        eos = equations.equation_of_state

        V, v1, v2, T = cons2thermo(u, equations)
    
        rho = u[1]
        p = pressure(V, T, eos)
        
        # critical compressibility
        Zc, Vc, Pc = critical_compressibility(eos)
        pc = Trixi.pressure(Vc, eos.Tc, eos)
        u0 = 25 # from initial condition

        return SVector(rho * Vc, v1 / u0, v2, p / pc)
    end

    Trixi.varnames(::typeof(cons2relative), ::NonIdealCompressibleEulerEquations2D) = 
        ("rho / rho_c", "u/u0", "v2", "p/pc")

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

    @show extrema(indicator.cache.alpha)

    using Plots
    pd = PlotData2D(sol.u[end], semi; solution_variables = cons2relative)
    p1 = plot(pd["rho / rho_c"], clims=(0.5, 2.5), dpi=400)
    if volume_integral == volume_integral_apec
        savefig(p1, "figs/apec_polydeg$(Trixi.polydeg(basis))_transcritical_mixing_rho_rho_c.png")
    else
        savefig(p1, "figs/ec_polydeg$(Trixi.polydeg(basis))_transcritical_mixing_rho_rho_c.png")
        # savefig(p1, "figs/ec_polydeg$(Trixi.polydeg(basis))_transcritical_mixing_rho_rho_c_long.png")
    end

    p2 = plot(ScalarPlotData2D(gamma, sol.u[end], semi), dpi=400)
    if volume_integral == volume_integral_apec
        savefig(p2, "figs/apec_polydeg$(Trixi.polydeg(basis))_transcritical_mixing_gamma.png")
    else
        savefig(p2, "figs/ec_polydeg$(Trixi.polydeg(basis))_transcritical_mixing_gamma.png")
        # savefig(p2, "figs/ec_polydeg$(Trixi.polydeg(basis))_transcritical_mixing_gamma_long.png")
    end

end # save_plots