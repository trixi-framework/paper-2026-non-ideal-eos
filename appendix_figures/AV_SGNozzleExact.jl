using StaticArrays
using Clapeyron, ForwardDiff
using StartUpDG
using OrdinaryDiffEq
using RecursiveArrayTools
using StaticArrays
using Trixi
using Plots
using MAT

include("SG_utils.jl")
function area(x)
    return 1 + 0.5 * cos(2*pi*x)
end

############## ############## ############## ############## 
############## Subsonic ############## ############## 
# subsonic case all around, location of throat x = 0.5
model = MySG();
equations = QuasiEuler1D(model)
gamma = model.gamma
cv = model.cv
Pinf = model.pInf
Pb = 0.5e6
P0 = 1e6;
T0 = 453.0;

PReff = (Pb + Pinf)/(P0 + Pinf)
As_Astar = area(1.0) / area(0.5)

f(M) = 1/(As_Astar * M) * ((1 + ((gamma - 1)/2) * 
    M^2)/((gamma + 1)/2))^(gamma +1)/(2 * (gamma - 1)) - 1
dfdM(M) = ForwardDiff.derivative(m -> f(m), M)

M = 1.0;
while (abs(f(M)) > 100*eps())
    M = M - f(M)/dfdM(M);
end

# solve for Mb (mach number at outlet)
s(M) = (1 + (gamma - 1)/2 *M^2)^(-gamma/(gamma - 1)) - PReff
dsdM(M) = ForwardDiff.derivative(m -> s(m), M)
M = 0.1
while (abs(s(M)) > 100*eps())
    M = M - s(M)/(dsdM(M))
end

cpr1 = (1 + (gamma - 1)/2 * M^2)^(-gamma/(gamma + 1))
Mb = M
Tb =  (T0)/(1 + (gamma - 1)/2 * Mb^2)
rho_b = density_pT(Pb, Tb, equations)
c_b = sqrt(gamma * (Pb + Pinf)/rho_b)
u_b = Mb * c_b
massflow = rho_b * u_b * area(1.0)

# knowing mass flow rate, by mass conservation we can find the others 
T(M) = T0/(1 + (gamma - 1)/2 *M^2)
P(M) = (P0 + Pinf)* (1 + (gamma - 1)/2 * M^2)^(-(gamma)/(gamma - 1)) - Pinf
rho(M) = (P(M) + Pinf)/((gamma - 1) * cv * T(M))
u_vel(M) = M * sqrt(gamma * (P(M) + Pinf)/rho(M))

# need pointwise root for each x, 
N = 3
M = 100
rd = RefElemData(Line(), SBP(), N)
(VX, ), EToV = uniform_mesh(Line(), M)

function domain_change(x) 
    a = 0
    b = 1.0;
    return (b - a)/2 * (x + 1) +a;
end

VX = domain_change.(VX)
md = MeshData((VX,), EToV, rd, is_periodic = false)
M_x = similar(md.x)
for e in axes(md.x, 2)
    for i in axes(md.x, 1)
        M = 0.06
        foundM(M) = area(md.x[i, e]) * rho(M) * u_vel(M) - massflow
        dDm(M) = ForwardDiff.derivative(m -> foundM(m), M)
        while (abs(foundM(M)) > 100000*eps())
            M = M - foundM(M)/ dDm(M)
            #@show abs(foundM(M))
        end
        M_x[i, e] = M
    end
end

## with Mach number, get the other variables
P_x = similar(md.x)
T_x = similar(md.x)
u_x = similar(md.x)
rho_x = similar(md.x)
for e in axes(md.x, 2)
    for i in axes(md.x, 1)
        T_x[i, e] = T0/(1 + (gamma - 1)/2 * M_x[i, e]^2);
        P_x[i, e] = (P0 + Pinf) /((1 + (gamma - 1)/2 * M_x[i, e]^2)^(gamma/(gamma - 1))) - Pinf
        rho_x[i, e] = (P_x[i, e] + Pinf)/((gamma - 1) * cv * T_x[i, e])
        u_x[i, e] = M_x[i, e] * sqrt(gamma * (P_x[i, e] + Pinf)/rho_x[i, e])
    end
end

N = size(md.x, 1)
x = [md.x[:, 1] ; vec(md.x[2:N, :])]
rho_data = [rho_x[:, 1] ; vec(rho_x[2:N, :])]
u_data = [u_x[:, 1] ; vec(u_x[2:N, :])]
p_data = [P_x[:, 1] ; vec(P_x[2:N, :])]
M_data = [M_x[:, 1] ; vec(M_x[2:N, :])]
T_data = [T_x[:, 1] ; vec(T_x[2:N, :])]
plot(x, p_data)
## write exact solution to .mat file. using MAT
data = Dict(
    "x" => x,
    "rho" => rho_data,
    "u" => u_data,
    "p" => p_data,
    "Mach" => M_data,
    "Temp" => T_data
)

# Write the dictionary to a .mat file
matwrite("exactWater.mat", data)

############## ############## ############## ############## 
############## Supersonic ############## ############## 
model = MySG(0.0, 2030e3, 1.0, 1.43, 1040.0)
equations = QuasiEuler1D(model)
gamma = model.gamma
cv = model.cv
Pinf = model.pInf
Pb = 0.5e6
P0 = 1e6;
T0 = 453.0;

# first guess shock point xs
function shockpoint(xs) 
    # preshock state (1)
    A_Ac = area(xs) / area(0.5)
    fM1(M1) = 1/(M1) * (2 / (gamma + 1) *(1 + (gamma - 1)/2 *M1^2))^
        ((gamma + 1)/(2 * (gamma - 1))) - A_Ac
    dfdM1(M1) = ForwardDiff.derivative(m -> fM1(m), M1)

    M1 = 1.1
    while (abs(fM1(M1)) > 100*eps()) 
        M1 = M1 - fM1(M1)/dfdM1(M1)
    end

    T1 = T0 /(1 + (gamma - 1)/ 2 * M1^2)
    P1 = (P0 + Pinf) * (1 + (gamma - 1)/2 * M1^2)^(-gamma/(gamma - 1)) - Pinf 
    rho1 = (P1 + Pinf)/((gamma - 1)*cv * T1)
    c1 = sqrt(gamma * (P1 + Pinf)/rho1)
    u1 = M1 * c1

    # now compute post shock state (2)
    M2 = sqrt((1 + (gamma - 1)/2 * M1^2)/(gamma * M1^2 - (gamma - 1)/2))
    P2 = ((2*gamma/(gamma + 1) * M1^2) - (gamma - 1)/(gamma + 1)) * (P1 + Pinf) - Pinf 
    T2 = ((2 * gamma)/(gamma + 1) * M1^2 - (gamma - 1)/(gamma + 1)) * 
    (1 + (gamma - 1)/2 * M1^2) / ((gamma + 1)^2/2 * M1^2) * T1
    rho2 = (P2 + Pinf)/((gamma - 1) * cv * T2)
    c2 = sqrt(gamma * (P2 + Pinf)/rho2)
    u2 = M2 * c2

    # new stagnation state after shock 
    T02 = T2 * (1 + (gamma - 1)/2 * M2^2)
    P02 = (P2 + Pinf) * (1 + (gamma - 1)/2 * M2^2)^(gamma/(gamma - 1))

    # compute exit mach, pressure, etc.
    A1_As = area(1)/area(xs)
    fMe(Me) = M2/Me * ((1 + (gamma - 1)/2 * Me^2)/(1 + (gamma - 1)/2 * M2^2))^
                ((gamma + 1)/(2 * (gamma - 1))) - A1_As
    dfdMe(Me) = ForwardDiff.derivative(m -> fMe(m), Me)

    Me = 0.5
    while (abs(fMe(Me)) > 100*eps()) 
        Me = Me - fMe(Me)/ dfdMe(Me)
    end
    Te = (T02/(1 + (gamma - 1)/2 * Me^2))
    Pe = (P02 + Pinf) * (1 + (gamma - 1)/2 * Me^2)^(-gamma/(gamma - 1)) - Pinf 

    # Pe must match Pb,
    return Pe - Pb
end

# first guess shock point xs
function shockpoint_helper(xs) 
    # preshock state (1)
    A_Ac = area(xs) / area(0.5)
    fM1(M1) = 1/(M1) * (2 / (gamma + 1) *(1 + (gamma - 1)/2 *M1^2))^
        ((gamma + 1)/(2 * (gamma - 1))) - A_Ac
    dfdM1(M1) = ForwardDiff.derivative(m -> fM1(m), M1)

    M1 = 1.1
    while (abs(fM1(M1)) > 100*eps()) 
        M1 = M1 - fM1(M1)/dfdM1(M1)
    end

    T1 = T0 /(1 + (gamma - 1)/ 2 * M1^2)
    P1 = (P0 + Pinf) * (1 + (gamma - 1)/2 * M1^2)^(-gamma/(gamma - 1)) - Pinf 
    rho1 = (P1 + Pinf)/((gamma - 1)*cv * T1)
    c1 = sqrt(gamma * (P1 + Pinf)/rho1)
    u1 = M1 * c1

    # now compute post shock state (2)
    M2 = sqrt((1 + (gamma - 1)/2 * M1^2)/(gamma * M1^2 - (gamma - 1)/2))
    P2 = ((2*gamma/(gamma + 1) * M1^2) - (gamma - 1)/(gamma + 1)) * (P1 + Pinf) - Pinf 
    rho2 = rho1 * ((gamma + 1)*M1^2)/(2 + (gamma - 1) * M1^2)
    T2 = (P2 + Pinf)/((gamma - 1) * cv * rho2)
    c2 = sqrt(gamma * (P2 + Pinf)/rho2)
    u2 = M2 * c2

    # new stagnation state after shock 
    T02 = T2 * (1 + (gamma - 1)/2 * M2^2)
    P02 = (P2 + Pinf) * (1 + (gamma - 1)/2 * M2^2)^(gamma/(gamma - 1))

    # compute exit mach, pressure, etc.
    A1_As = area(1)/area(xs)
    fMe(Me) = M2/Me * ((1 + (gamma - 1)/2 * Me^2)/(1 + (gamma - 1)/2 * M2^2))^
                ((gamma + 1)/(2 * (gamma - 1))) - A1_As
    dfdMe(Me) = ForwardDiff.derivative(m -> fMe(m), Me)

    Me = 0.5
    while (abs(fMe(Me)) > 100*eps()) 
        Me = Me - fMe(Me)/ dfdMe(Me)
    end
    return T02, P02, M2
end

dshockdxs(xs) = ForwardDiff.derivative(x -> shockpoint(x), xs)
xs = 0.7
iter = 0;
while (abs(shockpoint(xs))/Pb > 100* eps() && iter < 100)
    xs = xs - shockpoint(xs) / dshockdxs(xs)
    iter += 1
end

T02, P02, M2 = shockpoint_helper(xs)

# after getting shock, now get the solutioon profile
# upstream 0 < x < xc (throat point)
T01 = T0
P01 = P0
for e in axes(md.x, 2)
    for i in axes(md.x, 1)
        if (md.x[i, e] - 0.5 < 1000*eps())
            M = 0.2
            Astar = area(0.5)
            T0 = T01
            P0 = P01
        elseif (abs(md.x[i, e] - 0.5) <= 1000*eps())
            M = 1.0;
            T0 = T01
            P0 = P01
            Astar = area(0.5)
        elseif ((0.5 + 1000*eps() < md.x[i, e]) && (md.x[i, e] < xs - 1e-8))
            M = 2.0;
            T0 = T01
            P0 = P01
            Astar = area(0.5)
        elseif (abs(md.x[i, e] - xs) < 100*eps())
            M = M2
            Astar = area(0.5)
            T0 = T02;
            P0 = P02;
            continue;
        else
            Astar = area(xs) / (1/M2 * (2/(gamma+1) * (1 + (gamma-1)/2*M2^2))^((gamma+1)/(2*(gamma-1))))
            M = 0.2
            T0 = T02;
            P0 = P02;
        end
        r(M) = 1/M*(2/(gamma + 1) * (1 + (gamma - 1)/2 * M^2))^((gamma + 1)/(2 * (gamma - 1)))-
            area(md.x[i, e]) / Astar;
        drdM(M) = ForwardDiff.derivative(m -> r(m), M)
        while (abs(r(M)) > 100*eps()) 
            M = M - r(M) / drdM(M)
        end
        M_x[i, e] = M
        T_x[i, e] = T0/(1 + (gamma - 1)/2 * M_x[i, e]^2)
        P_x[i, e] = (P0 + Pinf) * (1 + (gamma - 1)/2 * M_x[i, e]^2)^(-gamma/(gamma - 1)) - Pinf 
        rho_x[i, e] = (P_x[i, e] + Pinf)/((gamma - 1) * cv * T_x[i, e])
        u_x[i, e] = M_x[i, e] * sqrt(gamma * (P_x[i, e] + Pinf)/rho_x[i, e])
    end
end

N = size(md.x, 1)
x = [md.x[:, 1] ; vec(md.x[2:N, :])]
rho_data = [rho_x[:, 1] ; vec(rho_x[2:N, :])]
u_data = [u_x[:, 1] ; vec(u_x[2:N, :])]
p_data = [P_x[:, 1] ; vec(P_x[2:N, :])]
M_data = [M_x[:, 1] ; vec(M_x[2:N, :])]
T_data = [T_x[:, 1] ; vec(T_x[2:N, :])]
## write exact solution to .mat file. using MAT
data = Dict(
    "x" => x,
    "rho" => rho_data,
    "u" => u_data,
    "p" => p_data,
    "Mach" => M_data,
    "Temp" => T_data
)

matwrite("exactSteam.mat", data)