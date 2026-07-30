# Appendix figure reproduction (stiffened gas)

Scripts in this folder reproduce **Figures 10–11** and the stiffened-gas parameter table in the appendix of the compiled paper (`paper/main.pdf` in the parent research repository). Run all commands from this directory.

## Environment setup

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Julia **1.12.6** was used to resolve [`Manifest.toml`](Manifest.toml).

## Figure and table map

| Paper reference | Script(s) | Notes |
|-----------------|-----------|-------|
| **SG parameter table** (appendix) | [`SG_utils.jl`](SG_utils.jl) (`MySG` defaults) | Table values are authored in the manuscript; constants match the code |
| **Figure 10** — liquid-water quasi-nozzle | [`AV_SGNozzleExact.jl`](AV_SGNozzleExact.jl), [`AV_SGQuasiNozzle.jl`](AV_SGQuasiNozzle.jl) | Water at `N=3, M=100` and `N=7, M=50` |
| **Figure 11** — steam quasi-nozzle | same scripts | Steam at `N=3, M=100` and `N=7, M=50` |

## Files

| File | Role |
|------|------|
| [`SG_utils.jl`](SG_utils.jl) | Stiffened-gas EOS, fluxes, and quasi-1D helpers |
| [`AV_SGNozzleExact.jl`](AV_SGNozzleExact.jl) | Analytical nozzle solution; writes `exactWater.mat` and `exactSteam.mat` |
| [`AV_SGQuasiNozzle.jl`](AV_SGQuasiNozzle.jl) | DG simulation and comparison plots against analytical/reference data |
| [`run_generate_figs.jl`](run_generate_figs.jl) | Driver: exact subprocess, then four `trixi_include` quasi cases |
| [`exactWater.mat`](exactWater.mat) | Committed analytical/reference profile for water |
| [`exactSteam.mat`](exactSteam.mat) | Committed analytical/reference profile for steam |

## Run commands

Exact references plus all four water/steam quasi-nozzle panel resolutions:

```bash
julia --project=. --threads=auto run_generate_figs.jl
```

This runs `AV_SGNozzleExact.jl` in a subprocess, then uses `trixi_include` to override `fluid`, `N`, and `M` in `AV_SGQuasiNozzle.jl` for:

| fluid | N | M   |
|-------|---|-----|
| Water | 3 | 100 |
| Water | 7 | 50  |
| Steam | 3 | 100 |
| Steam | 7 | 50  |

Or individually:

```bash
# Regenerate analytical .mat references (optional; committed copies are provided)
julia --project=. AV_SGNozzleExact.jl

# Single solve (defaults: water, N=3, M=100); override via trixi_include as in run_generate_figs.jl
julia --project=. AV_SGQuasiNozzle.jl
```

Outputs under `figs/` follow `pressureQuasi{Water|Steam}N{N}M{M}.png` and `MachQuasi{Water|Steam}N{N}M{M}.png` (plus density/entropy diagnostics).
