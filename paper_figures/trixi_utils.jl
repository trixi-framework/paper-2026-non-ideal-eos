using Trixi
using Trixi: StaticArrays
using .StaticArrays: @SMatrix
using LinearAlgebra: norm

function temperature_pV(p, V, eos::IdealGas)
    (; R) = eos
    T = p * V / R
    return T
end

function temperature_pV(p, V, eos::VanDerWaals)
    (; R, a, b) = eos
    T = (p + a / V^2) * (V - b) / R
    return T
end

function temperature_pV(p, V, eos::PengRobinson)
    T = eos.Tc
    tol = 100 * eps(typeof(T))
    dp = pressure(V, T, eos) - p
    iter = 1
    while abs(dp) / abs(p) > tol && iter < 100
        dp = pressure(V, T, eos) - p
        dpdT_V = ForwardDiff.derivative(T -> pressure(V, T, eos), T)
        T = max(tol, T - dp / dpdT_V)
        iter += 1
    end
    if iter == 100
        @warn "Solver for temperature(V, p) did not converge"
    end
    return T
end

function calc_barth_matrices(u_ll, u_rr, equations)

    # evaluate barth matrices at average state
    u = 0.5 * (u_ll + u_rr)

    eos = equations.equation_of_state
    rho, v1, p = cons2prim(u, equations)
    V, v1, T = cons2thermo(u, equations)

    # calc thermo states
    e = Trixi.energy_internal_specific(V, T, eos)
    H = e + p * V + 0.5 * v1^2 # total enthalpy

    # thermodynamic derivatives
    dpdT_V, dpdV_T = Trixi.calc_pressure_derivatives(V, T, eos)
    (; cv0, R, b) = eos

    # calculate ratio of specific heats
    K1 = Trixi.calc_K1(V, eos)
    d2aT = Trixi.peng_robinson_d2a(T, eos)
    cp0 = cv0 + R
    cv = cv0 - K1 * T * d2aT
    cp = cp0 - R - K1 * T * d2aT - T * dpdT_V^2 / dpdV_T
    # gamma = cp / cv

    kappa = dpdT_V / cv

    c = Trixi.speed_of_sound(V, T, eos)

    drhoe_drho_p = Trixi.drho_e_internal_drho_at_const_p(V, T, eos)

    # # trying different averaging techniques for `drhoe_drho_at_const_p` term
    # V_ll, _, T_ll = cons2thermo(u_ll, equations)
    # V_rr, _, T_rr = cons2thermo(u_rr, equations)
    # drhoe_drho_p_ll = Trixi.drho_e_internal_drho_at_const_p(V_ll, T_ll, eos)
    # drhoe_drho_p_rr = Trixi.drho_e_internal_drho_at_const_p(V_rr, T_rr, eos)
    # drhoe_drho_p = 2 * drhoe_drho_p_ll * drhoe_drho_p_rr / (drhoe_drho_p_ll + drhoe_drho_p_rr) # harmonic mean
    # drhoe_drho_p = 0.5 * (drhoe_drho_p_ll + drhoe_drho_p_rr)

    R = @SMatrix [1 1 1;
                  v1-c v1 v1+c;
                  H-v1 * c (0.5 * v1^2+drhoe_drho_p) H+v1 * c]

    # cp = cv - T * dpdT_V^2 / dpdV_T
    D = Diagonal(SVector(
                     2 * c^2 / (rho * T),
                     rho * c^4 / (kappa^2 * cp * T^2),
                     2 * c^2 / (rho * T)
                 ))
    #R̃ = R * Diagonal(inv.(sqrt.(diag(D))))

    if !(abs.(diag(D)) ≈ diag(D))
        @show diag(D)
        # @show rho, v1, T
    end

    # code hardening to avoid d_2 ≤ 0
    D_eps = SVector(D.diag[1], max(1e-6, D.diag[2]), D.diag[3])
    R̃ = R * Diagonal(inv.(sqrt.(D_eps)))

    lam = SVector(v1 - c, v1, v1 + c)
    Λplus = Diagonal(abs.(lam))
    return R̃, Λplus
end

@inline function dissipation_barth(u_ll, u_rr, orientation::Int, equations)
    v_ll = cons2entropy(u_ll, equations)
    v_rr = cons2entropy(u_rr, equations)
    dv = v_rr - v_ll
    R, Lambda = calc_barth_matrices(u_ll, u_rr, equations)
    # Lambda = max_abs_speed(u_ll, u_rr, 1, equations) # scalar dissipation
    return -0.5f0 * R * Lambda * R' * dv
end

@inline function dissipation_barth_with_sensor(u_ll, u_rr, orientation::Int, equations)
    v_ll = cons2entropy(u_ll, equations)
    v_rr = cons2entropy(u_rr, equations)
    dv = v_rr - v_ll
    R, Lambda = calc_barth_matrices(u_ll, u_rr, equations)

    # p_ll = pressure(u_ll, equations)
    # p_rr = pressure(u_rr, equations)
    # sensor_pressure = min(1, 2 * abs(p_ll - p_rr) / min(p_ll, p_rr) + 0.1)

    # make sensor proportional to [ρ] with some small fallback
    rho_ll = Trixi.density(u_ll, equations)
    rho_rr = Trixi.density(u_rr, equations)
    sensor_pressure = min(1, abs(rho_ll - rho_rr) / min(rho_ll, rho_rr) + 0.1)

    # scale the contact wave eigenvalue
    sensor = SVector(1, sensor_pressure, 1)
    scaled_Lambda = Diagonal(sensor .* Lambda.diag)
    return -0.5 * R * scaled_Lambda * R' * dv
end


@inline function drho_e_dp_at_const_rho(V, T, eos::Trixi.AbstractEquationOfState)
    rho = inv(V)
    dpdT_V, _dpdV_T = Trixi.calc_pressure_derivatives(V, T, eos)
    c_v = Trixi.heat_capacity_constant_volume(V, T, eos)

    # (∂(ρe)/∂p)|ρ = ρ c_v / (∂p/∂T)|V
    return (rho * c_v) / dpdT_V
end


@inline function flux_epec_ranocha(u_ll, u_rr, orientation::Int, equations::NonIdealCompressibleEulerEquations1D)
    eos = equations.equation_of_state
    V_ll, v1_ll, T_ll = cons2thermo(u_ll, equations)
    V_rr, v1_rr, T_rr = cons2thermo(u_rr, equations)

    rho_ll = u_ll[1]
    rho_rr = u_rr[1]
    e_internal_ll = Trixi.energy_internal_specific(V_ll, T_ll, eos)
    e_internal_rr = Trixi.energy_internal_specific(V_rr, T_rr, eos)
    rho_e_internal_ll = rho_ll * e_internal_ll
    rho_e_internal_rr = rho_rr * e_internal_rr
    p_ll = pressure(V_ll, T_ll, eos)
    p_rr = pressure(V_rr, T_rr, eos)

    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    rho_e_internal_avg = 0.5f0 * (rho_e_internal_ll + rho_e_internal_rr)
    p_v1_avg = 0.5f0 * (p_ll * v1_rr + p_rr * v1_ll)

    # chain rule from Terashima
    drho_e_internal_drho_p_ll = Trixi.drho_e_internal_drho_at_const_p(V_ll, T_ll, eos)
    drho_e_internal_drho_p_rr = Trixi.drho_e_internal_drho_at_const_p(V_rr, T_rr, eos)
    drho_e_internal_drho_p_avg = 0.5f0 * (drho_e_internal_drho_p_ll + drho_e_internal_drho_p_rr)
    drho_e_internal_drho_p_rho_avg = 0.5f0 * (drho_e_internal_drho_p_ll * rho_ll + drho_e_internal_drho_p_rr * rho_rr)

    rho_e_jump = rho_e_internal_rr - rho_e_internal_ll
    rho_jump = rho_rr - rho_ll
    p_jump = p_rr - p_ll
    drho_e_internal_drho_p_jump = drho_e_internal_drho_p_rr - drho_e_internal_drho_p_ll
    drho_e_dp_at_const_rho_ll = drho_e_dp_at_const_rho(V_ll, T_ll, eos)
    drho_e_dp_at_const_rho_rr = drho_e_dp_at_const_rho(V_rr, T_rr, eos)
    drho_e_dp_at_const_rho_avg = 0.5f0 * (drho_e_dp_at_const_rho_ll +
                                          drho_e_dp_at_const_rho_rr)
    num = (rho_e_jump - drho_e_internal_drho_p_avg * rho_jump - drho_e_dp_at_const_rho_avg * p_jump)
    den = drho_e_internal_drho_p_jump
    rho_avg = rho_avg - num * den / (den^2 + eps(typeof(den)))

    rho_e_v1_avg = (rho_e_internal_avg + drho_e_internal_drho_p_avg * rho_avg - drho_e_internal_drho_p_rho_avg) *
                   v1_avg

    # check Ranocha condition
    check_condition = false
    if check_condition == true
        lhs = (drho_e_internal_drho_p_rr - drho_e_internal_drho_p_ll) * rho_avg
        rhs = (drho_e_internal_drho_p_rr * rho_rr - drho_e_internal_drho_p_ll * rho_ll) -
            (rho_e_internal_rr - rho_e_internal_ll)
        @show lhs - rhs
    end


    # Ignore orientation since it is always "1" in 1D
    f_rho = rho_avg * v1_avg
    f_rho_v1 = rho_avg * v1_avg * v1_avg + p_avg
    f_rho_E = rho_e_v1_avg + rho_avg * 0.5f0 * (v1_ll * v1_rr) * v1_avg + p_v1_avg

    return SVector(f_rho, f_rho_v1, f_rho_E)
end

@inline function flux_coppola(u_ll, u_rr, orientation::Int,
                              equations::NonIdealCompressibleEulerEquations1D)
    eos = equations.equation_of_state
    V_ll, v1_ll, T_ll = cons2thermo(u_ll, equations)
    V_rr, v1_rr, T_rr = cons2thermo(u_rr, equations)

    rho_ll = u_ll[1]
    rho_rr = u_rr[1]
    e_internal_ll = Trixi.energy_internal_specific(V_ll, T_ll, eos)
    e_internal_rr = Trixi.energy_internal_specific(V_rr, T_rr, eos)
    rho_e_internal_ll = rho_ll * e_internal_ll
    rho_e_internal_rr = rho_rr * e_internal_rr
    p_ll = pressure(V_ll, T_ll, eos)
    p_rr = pressure(V_rr, T_rr, eos)

    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    p_v1_avg = 0.5f0 * (p_ll * v1_rr + p_rr * v1_ll)

    # Ignore orientation since it is always "1" in 1D

    # coppola's new flux
    drho_e_internal_drho_p_ll = Trixi.drho_e_internal_drho_at_const_p(V_ll, T_ll, eos)
    drho_e_internal_drho_p_rr = Trixi.drho_e_internal_drho_at_const_p(V_rr, T_rr, eos)
    drho_e_internal_drho_p_avg = 0.5f0 * (drho_e_internal_drho_p_ll + drho_e_internal_drho_p_rr)

    # (drho_e/drho)_p = rho * (de/drho)_p + e
    # --> rho^2 * (de/drho)_p = rho * (drho_e/drho)_p - rho_e
    rho_squared_de_drho_p_ll = drho_e_internal_drho_p_ll * rho_ll - rho_e_internal_ll
    rho_squared_de_drho_p_rr = drho_e_internal_drho_p_rr * rho_rr - rho_e_internal_rr
    rho_squared_de_drho_p_avg = 0.5f0 * (rho_squared_de_drho_p_ll + rho_squared_de_drho_p_rr)

    use_EPEP = true
    if use_EPEP == true
        # Exactly PEP: use rho_avg = [rho^2 * (de/drho)_p] / [drho_e_internal_drho_p]
        num = rho_squared_de_drho_p_rr - rho_squared_de_drho_p_ll
        denom = drho_e_internal_drho_p_rr - drho_e_internal_drho_p_ll
        if abs(denom) < 100 * eps(typeof(denom))
            rho_lambda_avg = rho_avg
        else
            rho_lambda_avg = num / denom
        end
        # rho_lambda_avg = Trixi.regularized_ratio(num, denom)


        # drho_e_internal_drho_p_avg * f_rho - v1_avg * 0.5f0 * (rho_squared_de_drho_p_ll + rho_squared_de_drho_p_rr)
        f_rho = rho_lambda_avg * v1_avg
        rho_e_internal_corrected =
            (drho_e_internal_drho_p_avg * rho_lambda_avg - rho_squared_de_drho_p_avg)
    else
        f_rho = rho_avg * v1_avg
        rho_e_internal_corrected =
            (drho_e_internal_drho_p_avg * rho_avg - rho_squared_de_drho_p_avg)
        rho_lambda_avg = rho_avg # for checking condition
    end

    # check Ranocha condition
    check_condition = false
    if check_condition == true
        # @show (drho_e_internal_drho_p_rr - drho_e_internal_drho_p_ll) * rho_lambda_avg
        # @show (drho_e_internal_drho_p_rr * rho_rr - drho_e_internal_drho_p_ll * rho_ll) -
        #     (rho_e_internal_rr - rho_e_internal_ll)
        lhs = (drho_e_internal_drho_p_rr - drho_e_internal_drho_p_ll) * rho_lambda_avg
        rhs = (drho_e_internal_drho_p_rr * rho_rr - drho_e_internal_drho_p_ll * rho_ll) -
            (rho_e_internal_rr - rho_e_internal_ll)
        @show lhs - rhs
    end

    f_rho_v1 = f_rho * v1_avg + p_avg
    f_rho_E = rho_e_internal_corrected * v1_avg + f_rho * 0.5f0 * (v1_ll * v1_rr) + p_v1_avg

    return SVector(f_rho, f_rho_v1, f_rho_E)
end

# returns both u and x
function sol_and_coordinates(sol, index=length(sol.u))
    semi = sol.prob.p
    x = semi.cache.elements.node_coordinates
    x = reshape(x, size(x, 2), size(x, 3))
    u = sol.u[index]
    u = wrap_and_reshape(u, semi)
    return u, x
end

function wrap_and_reshape(u, semi)
    u = Trixi.wrap_array_native(u, semi)
    u = reinterpret(SVector{3, Float64}, u)
    return reshape(u, size(u, 2), size(u, 3)) # assume first dim = 1
end

function calc_solution_norm_1d(u, semi)
    u = wrap_and_reshape(u, semi)
    J = inv.(semi.cache.elements.inverse_jacobian)
    w = semi.solver.basis.weights
    wJ = w * J'
    return sqrt(sum(wJ .* norm.(u).^2))
end

function calc_error_1d(sol; normalize=true)
    error = map(sol.u, sol.t) do u, t
        x = semi.cache.elements.node_coordinates
        x = reshape(x, size(x, 2), size(x, 3))
        u = wrap_and_reshape(u, semi)
        J = inv.(semi.cache.elements.inverse_jacobian)
        w = semi.solver.basis.weights
        wJ = w * J'

        solution_norm = sqrt(sum(wJ .* norm.(initial_condition.(x, t, equations)).^2))
        if normalize == true
            return sqrt(sum(wJ .* norm.(u - initial_condition.(x, t, equations)).^2)) / solution_norm
        else
            return sqrt(sum(wJ .* norm.(u - initial_condition.(x, t, equations)).^2))
        end
    end
    return error
end

# here, field is a function of the form field(u, equations) -> scalar
function calc_field_error_1d(sol; field=pressure)
    error = map(sol.u, sol.t) do u, t
        x = semi.cache.elements.node_coordinates
        x = reshape(x, size(x, 2), size(x, 3))
        u = wrap_and_reshape(u, semi)
        J = inv.(semi.cache.elements.inverse_jacobian)
        w = semi.solver.basis.weights
        wJ = w * J'

        solution_norm = sqrt(sum(wJ .* field.(initial_condition.(x, t, equations), equations).^2))
        return sqrt(sum(wJ .* (field.(u, equations) - field.(initial_condition.(x, t, equations), equations)).^2)) / solution_norm
    end
    return error
end

function calc_entropy_1d(sol)
    entropy = map(sol.u, sol.t) do u, t
        u = wrap_and_reshape(u, semi)
        J = inv.(semi.cache.elements.inverse_jacobian)
        w = semi.solver.basis.weights
        wJ = w * J'
        return sum(wJ .* Trixi.entropy.(u, equations))
    end
    return (entropy .- entropy[1]) / maximum(abs.(entropy))
end