# Run appendix figure scripts: exact SG nozzle references, then all
# water/steam quasi-nozzle panel resolutions via trixi_include.
#
# Usage (from appendix_figures/):
#   julia --project=. --threads=auto run_generate_figs.jl

cd(@__DIR__)

julia = Base.julia_cmd()
project = @__DIR__

println("="^64)
println("Running AV_SGNozzleExact.jl")
println("="^64)
run(`$julia --project=$project --threads=auto AV_SGNozzleExact.jl`)

using Trixi

const quasi_cases = (
    (fluid = "Water", N = 3, M = 100),
    (fluid = "Water", N = 7, M = 50),
    (fluid = "Steam", N = 3, M = 100),
    (fluid = "Steam", N = 7, M = 50),
)

for case in quasi_cases
    println("="^64)
    println("Running AV_SGQuasiNozzle.jl  fluid=$(case.fluid)  N=$(case.N)  M=$(case.M)")
    println("="^64)
    trixi_include("AV_SGQuasiNozzle.jl";
                  fluid = case.fluid,
                  N = case.N,
                  M = case.M)
end

println("Appendix figure scripts finished. See figs/ and exact*.mat outputs.")
