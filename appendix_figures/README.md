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
| **Figure 10** — liquid-water quasi-nozzle | [`AV_SGNozzleExact.jl`](AV_SGNozzleExact.jl), [`AV_SGQuasiNozzle.jl`](AV_SGQuasiNozzle.jl) | Pressure and Mach panels for water |
| **Figure 11** — steam quasi-nozzle | same scripts with steam model toggles | Two rows: `N=3, M=100` and `N=7, M=50` |

## Files

| File | Role |
|------|------|
| [`SG_utils.jl`](SG_utils.jl) | Stiffened-gas EOS, fluxes, and quasi-1D helpers |
| [`AV_SGNozzleExact.jl`](AV_SGNozzleExact.jl) | Analytical nozzle solution; writes `exactWater.mat` and `exactSteam.mat` |
| [`AV_SGQuasiNozzle.jl`](AV_SGQuasiNozzle.jl) | DG simulation and comparison plots against analytical/reference data |
| [`exactWater.mat`](exactWater.mat) | Committed analytical/reference profile for water |
| [`exactSteam.mat`](exactSteam.mat) | Committed analytical/reference profile for steam |

## Run commands

Exact reference plus the default quasi-nozzle solve (water, `N=3`, `M=100`) can be run together:

```bash
julia --project=. --threads=auto run_generate_figs.jl
```

Or individually:

```bash
# Regenerate analytical .mat references (optional; committed copies are provided)
julia --project=. AV_SGNozzleExact.jl

# Solve and write publication PNGs (default: water, N=3, M=100)
julia --project=. AV_SGQuasiNozzle.jl
```

Default outputs under `figs/` (with the water/`N=3`/`M=100` settings):

- `figs/pressureQuasiWaterN3M100.png`
- `figs/MachQuasiWaterN3M100.png`

Pressure and Mach filenames match the manuscript naming convention (`pressureQuasi{Water|Steam}N{N}M{M}.png`, `MachQuasi{Water|Steam}N{N}M{M}.png`). Density and entropy PNGs are also written for diagnostics.

## Manual configuration for Figure 10–11

[`AV_SGQuasiNozzle.jl`](AV_SGQuasiNozzle.jl) ships with the **water** model, `N = 3`, and `M = 100`. Edit before running:

1. Comment/uncomment the water vs steam `model` constructor.
2. Set `fluid = "Water"` or `fluid = "Steam"` to match (controls output filenames).
3. For water Figure 10, set `N = 7`, `M = 50` (paper panel).
4. For steam Figure 11, run twice: (`N=3`, `M=100`) and (`N=7`, `M=50`).
5. For steam, also point `matread(...)` at `exactSteam.mat` instead of `exactWater.mat`.

## Remaining limitations

Water vs steam, polynomial degree, and mesh resolution still require manual edits rather than command-line arguments.
