# Run all generate_*.jl scripts (Figures 1–5) from this directory.
#
# Usage (from paper_figures/):
#   julia --project=. --threads=auto run_generate_figs.jl

const scripts = (
    "generate_EPEC_vdw_figs.jl",
    "generate_APEC_LxF_pr_figs.jl",
    "generate_APEC_EC_pr_figs.jl",
    "generate_pr_shock_figs.jl",
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

println("All generate_*.jl scripts finished. Figures written under figs/.")
