# Minimal verification of Appendix B entropy-Hessian factorizations
# (∂u/∂q, ∂q/∂v) for Trixi van der Waals and Peng-Robinson EOS in 1D.
# The 2D case adds one velocity row/column with the same structure.
#
# Run: julia --project=. test_appendix_B_eos.jl

using ForwardDiff
using LinearAlgebra
using StaticArrays
using Trixi

"""
    dudq(q, equations)

Analytic ∂u/∂q for q = (V, v₁, T) from Appendix B (1D).
"""
function dudq(q, equations::NonIdealCompressibleEulerEquations1D)
    eos = equations.equation_of_state
    V, v1, T = q
    rho = inv(V)
    rho2 = rho^2
    e = energy_internal_specific(V, T, eos)
    dpdT_V, _ = Trixi.calc_pressure_derivatives(V, T, eos)
    p = pressure(V, T, eos)
    dedV_T = T * dpdT_V - p
    dedT_V = Trixi.heat_capacity_constant_volume(V, T, eos)
    drhoEdV = (e + 0.5 * v1^2) * (-rho2) + dedV_T * rho
    drhoEdT = rho * dedT_V
    # column-major SMatrix
    return SMatrix{3, 3}(-rho2, -v1 * rho2, drhoEdV,
                         0, rho, rho * v1,
                         0, 0, drhoEdT)
end

"""
    dqdv(q, equations)

Analytic ∂q/∂v for q = (V, v₁, T) from Appendix B (1D), with
a = (∂(g/T)/∂V)_T⁻¹ and b = ½(v₁²) - T² (∂(g/T)/∂T)_V.
"""
function dqdv(q, equations::NonIdealCompressibleEulerEquations1D)
    eos = equations.equation_of_state
    V, v1, T = q
    dpdT_V, dpdV_T = Trixi.calc_pressure_derivatives(V, T, eos)
    p = pressure(V, T, eos)
    e = energy_internal_specific(V, T, eos)
    h = e + p * V
    dginvTdV_T = (V / T) * dpdV_T
    dginvTdT_V = (V * T * dpdT_V - h) / T^2
    inva = inv(dginvTdV_T)
    b = dginvTdT_V + 0.5 * v1^2 / T^2
    # column-major SMatrix
    return SMatrix{3, 3}(inva, 0, 0,
                         v1 * inva, T, 0,
                         (v1^2 - b * T^2) * inva, v1 * T, T^2)
end

function check_appendix_B(eos, q; label)
    equations = NonIdealCompressibleEulerEquations1D(eos)
    rho = q[1]
    V = inv(rho)
    v1 = q[2]
    T = q[3]
    u = thermo2cons(q, equations)

    A0 = ForwardDiff.hessian(u_ -> entropy(u_, equations), u)
    A1 = ForwardDiff.jacobian(u_ -> flux(u_, 1, equations), u)
    du_dv_analytic = dudq(q, equations) * dqdv(q, equations)
    du_dv_fd = inv(A0)

    @assert cons2entropy(u, equations) ≈
            ForwardDiff.gradient(u_ -> entropy(u_, equations), u)

    @assert Matrix(dudq(q, equations)) ≈
            ForwardDiff.jacobian(q_ -> thermo2cons(q_, equations), q)

    @assert norm(Matrix(du_dv_analytic) - du_dv_fd) <
            100 * eps() * norm(du_dv_fd)

    @assert all(eigvals(Matrix(A0)) .> eps())

    @assert norm(A0 * A1 - (A0 * A1)') < 100 * eps() * norm(A0 * A1)

    println("PASS: Appendix B checks for ", label)
    return nothing
end

# Trixi N₂ defaults; fixed admissible supercritical / safe states (V > b, cᵥ > 0)
q = SVector(2.0, -0.5, 300.0)
check_appendix_B(VanDerWaals(), q; label="VanDerWaals")

q = SVector(1.7, -0.1, 300.0)
check_appendix_B(PengRobinson(), q; label="PengRobinson")

println("All Appendix B EOS checks passed.")
