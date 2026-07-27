const FIGDIR = joinpath(@__DIR__, "figs")
if !isdir(FIGDIR)
    mkpath(FIGDIR)
end

using Trixi
using MAT
using StaticArrays
using Clapeyron, ForwardDiff
using StartUpDG
using OrdinaryDiffEq
using OrdinaryDiffEqSSPRK
using LinearAlgebra
using RecursiveArrayTools
using StaticArrays
using Trixi
using Plots
using LaTeXStrings

# Match paper_figures generate_*_figs.jl font sizes
default(legendfontsize=12, guidefontsize=14, tickfontsize=12, titlefontsize=14, dpi=400)

include("SG_utils.jl")
regularized_ratio(a, b) = a * b / (1e-15 + b^2)

function area(x)
    return 1 + 0.5 * cos(2*pi*x)
end

p_inflow() = 1e6;
p_outflow() = 0.5e6;

function initial_condition_nozzle(x, t, equations::QuasiEuler1D) 

    T0 = 453.0 # Kelvin
    p = (p_outflow() - p_inflow())*x + p_inflow()
    rho_inflow = density_pT(p, T0, equations)
    u_inflow = 0.0

    A = area(x)
    rho_e = rho_e_rhoP(rho_inflow, p, equations)
    rho_E = rho_e + 0.5 * rho_inflow * u_inflow^2;
    return A* SVector(rho_inflow, rho_inflow * u_inflow, rho_E, 1)

end

# if subsonic, the specified outlet pressure 
# msut be the ambient pressure (p_star)
# for the nozzle problem
# https://www.sciencedirect.com/science/article/pii/S0029549310005157?ref=pdf_download&fr=RR-2&rr=9bcddc6449a3290b
function nozzle_outlet_subsonic(u_ll, p_star, equations::QuasiEuler1D)
    model = equations.eos 
    (; pInf, gamma, q) = model
    a_rho, a_rho_v, a_rho_E, a = u_ll
    v_ll = a_rho_v / a_rho 
    p_ll = pressure(u_ll, equations)
    c_ll = speed_of_sound(u_ll, equations) 
    rho_ll = a_rho / a

    rho_star = rho_ll + (p_star - p_ll)/c_ll^2
    z_ll = rho_ll * c_ll
    v_star = v_ll + (p_ll - p_star)/z_ll

    rho_e_star = (p_star + gamma * pInf)/(gamma - 1) + rho_star * q
    a_rho_E_star = a * (rho_e_star + 0.5 *rho_star * v_star^2)
    return SVector(rho_star * a, 
        a * rho_star * v_star, 
        a_rho_E_star, a)
end

# inlet BC stagnation. It seems like it's a isentropic and isenthalpic 
# process which yields various equalities to solve for the
# the surface variables. see above reference as well
# this is going to follow convention of the paper where
# u_rr is actually the surface value (uM) and we need to find u_star
# which is (uP) in DG notation
# 0 subscript represents x=-∞ stagnation values which are given to us
function inlet_stagnation_nozzle(u_rr, 
        p_0, T_0, equations)
    (; pInf, gamma, q, cv) = equations.eos
    rho_0 = density_pT(p_0, T_0, equations)
    a_rho, a_rho_v, a_rho_E, a = u_rr
    rho_rr = a_rho / a;
    c_rr = speed_of_sound(u_rr, equations)
    p_rr = pressure(u_rr, equations)
    vel_rr = a_rho_v/a_rho
    z_rr = rho_rr * c_rr
    g_ratio = (gamma/(gamma - 1))
    K = (p_0 + pInf)/(rho_0^gamma)
    H_bar = gamma *(p_0 + pInf)/(rho_0 *(gamma - 1))
    
    #following the reference convention
    # we use u_star as the velocity at surface (uP point)
    f(u_star) = K * ((1/g_ratio * (1/K))^(g_ratio)) *
        (H_bar- 0.5 * u_star^2)^g_ratio  - pInf - p_rr + z_rr * (vel_rr - u_star)
    
    dfdu(u_star) = -(K * ((1/g_ratio * (1/K))^(g_ratio)
            * (g_ratio) * u_star *
            (H_bar - 0.5 * u_star^2)^(g_ratio/gamma)) + z_rr)
    
    u_star = vel_rr
    iter = 0
    diff = f(u_star)
    while (abs(diff) > 100*eps() && iter < 100) 
        u_star = u_star - f(u_star)/(dfdu(u_star))
        diff = f(u_star);
        iter+= 1
    end
    p_star = p_rr - z_rr *(vel_rr - u_star)
    rho_star = ((p_star + pInf)/K)^(1/gamma)
    rho_e_star = (p_star + gamma * pInf)/(gamma - 1) +rho_star * q
    return SVector(a * rho_star, a * rho_star * u_star, 
        a * (rho_e_star + 0.5 * rho_star * u_star^2), a)
end

function rhs!(du_voa, u_voa, params, t)
    (; rd, md, invMQTr, equations, C12, Qr_skew) = params
    du = parent(du_voa)
    u = parent(u_voa)
    # calc volume residual
     
    fill!(du, zero(eltype(du)))
    for e in axes(u, 2)
        for i in axes(u, 1), j in axes(u, 1)
            du[i, e] += Qr_skew[i, j] * flux_nonsym(u[i, e], u[j, e], equations)
        end
      @. du[:, e] /= rd.M.diag
    end
    uM = rd.Vf * u
    uP = uM[md.mapP]

    uP[1] = inlet_stagnation_nozzle(uM[1], p_inflow(), 453.0, equations)
    uP[end] = nozzle_outlet_subsonic(uM[end], p_outflow(), equations)

    # calc DG gradient
    v = cons2entropy.(u, equations)

    vM = cons2entropy.(uM, equations)
    vP = cons2entropy.(uP, equations)
    interface_flux = @. 0.5 * (vP - vM) * (1 - C12 * sign(md.nxJ)) * md.nxJ
    theta = (md.rxJ .* (rd.Dr * v) + rd.LIFT * interface_flux) ./ md.J

    # calculate AV coefficient
    sigma = similar(theta)
    for i in eachindex(sigma)
        invA0 = inv(ForwardDiff.jacobian(u -> cons2entropy(u, equations), u[i])[1:3, 1:3])
        sigma[i] = [invA0 zeros(3, 1);
                            zeros(1, 3) 0] * theta[i]
    end
    dissipation = sum(rd.M * (md.J .* dot.(sigma, theta)), dims=1)
    delta = sum(dot.(v, rd.M * du), dims=1) + sum((@. psi(uM, equations) * md.nx), dims=1)
    epsilon = vec(@. regularized_ratio(-min(0, delta), dissipation))
    epsilon .= 0.0
    sigma *= Diagonal(epsilon)
    sigmaM = rd.Vf * sigma
    sigmaP = sigmaM[md.mapP]
    interface_flux = @. 0.5 * (sigmaP - sigmaM) * (1 + C12 * sign(md.nxJ)) * md.nxJ
    du .-= md.rxJ .* (rd.Dr * sigma) + rd.LIFT * interface_flux

    du .+= rd.LIFT * @. flux_LXF_SG(uM, uP, md.nx, equations)

    @. du = -du / md.J
    return du
end

function domain_change(x) 
    a = 0
    b = 1.0;
    return (b - a)/2 * (x + 1) +a;
end

N = parse(Int, get(ENV, "N", "3"))
M = parse(Int, get(ENV, "M", "100"))
fluid = get(ENV, "FLUID", "Water")
rd = RefElemData(Line(), SBP(), N)
(VX, ), EToV = uniform_mesh(Line(), M)

VX = domain_change.(VX)
md = MeshData((VX,), EToV, rd, is_periodic = false)
 
if fluid == "Water"
    model = MySG()
    exact_file = "exactWater.mat"
elseif fluid == "Steam"
    model = MySG(0.0, 2030e3, 1.0, 1.43, 1040.0)
    exact_file = "exactSteam.mat"
else
    error("Unknown FLUID=$(fluid); expected Water or Steam")
end

equations = QuasiEuler1D(model)

initial_condition = initial_condition_nozzle

(Qrh,), VhP, Ph = hybridized_SBP_operators(rd)
invMQTr = rd.M \ (rd.Dr' * rd.M)
Qr_skew = rd.M * rd.Dr - rd.Dr' * rd.M

volume_flux = flux_central
#volume_flux = flux_central_terashima_VT

params = (; rd, md, invMQTr, C12 = 1, initial_condition, equations, volume_flux,
            Qr_skew, Qrh_skew=Qrh - Qrh', VhP, Vh = VhP * rd.Vq, Ph,
            epsilon_save = Float64[], t_save = Float64[], unorm_save = Float64[])

rq, wq = gauss_quad(0, 0, rd.N);
Vq = vandermonde(Line(), rd.N, rq) / rd.VDM;
Pq = (Vq' * diagm(wq) * Vq) \ (Vq' * diagm(wq));
u = Pq * initial_condition.(Vq * md.x, 0.0, equations);
## plot initial

tspan = (0.0, 0.3)
ode = ODEProblem(rhs!, VectorOfArray(u), tspan, params)
println("Running ode solver")
sol = OrdinaryDiffEq.solve(ode, SSPRK43(), 
            abstol=1e-8, reltol=1e-6, 
            #adaptive=false, dt = 1e-7, #dt = 0.1 * estimate_h(rd, md) / (lambda * (2 * rd.N + 1)), 
            saveat=LinRange(tspan..., 100), 
            callback=AliveCallback(alive_interval=200))
u = parent(sol.u[end])

pad_nans(u) = vec([u; fill(NaN, 1, size(u, 2))])
cell_avg(u, rd) = vec(sum(Diagonal(rd.wq) * (rd.Vq * u), dims=1) ./ 2.0)

# load exact solution 
data = matread(exact_file)
d = size(data["x"])
idxs = 1:10:d[1]
idxs = vcat(1:10:200, 205:6:304)

mkpath("figs")
tag = "Quasi$(fluid)N$(N)M$(M)"

# plot density
plot(vec(md.x), vec(getindex.(u, 1) ./ getindex.(u, 4)), linewidth=4, leg=false)
# plot average density
plot!(vec(cell_avg(md.x, rd)), vec(cell_avg(getindex.(u, 1) ./ getindex.(u, 4), rd)),
    linewidth=2.5, leg=false)
if fluid == "Water"
    scatter!(data["x"][idxs], data["rho"][idxs],
        yticks=901.2:-0.2:899.6, ylims=(899.5, 901.2), ms=3, leg=false)
else
    scatter!(data["x"][idxs], data["rho"][idxs], ms=3,
        yticks=5.0:-0.5:0.55, ylims=(0.55, 5.0), leg=false)
end
xlabel!(L"x")
ylabel!("Density " * L"(kg/m^3)")
savefig(joinpath(FIGDIR, "density$(tag).png"))

# pressure plot 
plot(vec(md.x), vec(pressure.(u, equations)), linewidth=4, leg=false)
plot!(vec(cell_avg(md.x, rd)), vec(cell_avg(pressure.(u, equations), rd)),
    linewidth=2.5, leg=false)
scatter!(data["x"][idxs], data["p"][idxs], ms=3, leg=false)
xlabel!(L"x")
ylabel!("Pressure (Pa)")
savefig(joinpath(FIGDIR, "pressure$(tag).png"))

# Ma Plot
vel_final = getindex.(u, 2) ./ getindex.(u,1)
speedsound_final = speed_of_sound.(u, equations);
plot(vec(md.x), vec(vel_final ./ speedsound_final), linewidth=4, leg=false)
plot!(vec(cell_avg(md.x, rd)), vec(cell_avg(vel_final ./ speedsound_final, rd)),
    linewidth=2.5, leg=false)
scatter!(data["x"][idxs], data["Mach"][idxs], ms=3, leg=false)
xlabel!(L"x")
ylabel!("Ma")
savefig(joinpath(FIGDIR, "Mach$(tag).png"))


## calculate change in entropy 
total_entropy = zeros(size(sol.t))
for (i, t) in enumerate(sol.t)
    total_entropy[i] = sum(md.wJq .* entropy.(parent(sol.u[i]), equations))
end
plot(sol.t, total_entropy, legend=false)
xlabel!(L"t")
ylabel!("Entropy")
savefig(joinpath(FIGDIR, "entropy$(tag).png"))
