# Run appendix figure scripts (stiffened-gas exact + default quasi-nozzle) from this directory.
#
# Usage (from appendix_figures/):
#   julia --project=. --threads=auto run_generate_figs.jl
#
# Note: AV_SGQuasiNozzle.jl defaults to water, N=3, M=100. Other Figure 10–11
# panels still require the manual toggles documented in README.md.

const scripts = (
    "AV_SGNozzleExact.jl",
    "AV_SGQuasiNozzle.jl",
)

cd(@__DIR__)

julia = Base.julia_cmd()
project = @__DIR__

for script in scripts
    println("="^64)
    println("Running ", script)
    println("="^64)
    run(`$julia --project=$project --threads=auto $script`)
end

println("Appendix figure scripts finished. See figs/ and exact*.mat outputs.")
