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

```bash
# Regenerate analytical .mat references (optional; committed copies are provided)
julia --project=. AV_SGNozzleExact.jl

# Interactive Plots comparison (default: water, N=3, M=100)
julia --project=. AV_SGQuasiNozzle.jl
```

## Manual configuration for Figure 11 (steam)

[`AV_SGQuasiNozzle.jl`](AV_SGQuasiNozzle.jl) ships with the **water** model active. For steam panels, edit the script before running:

1. Comment the default `model = MySG()` water constructor.
2. Uncomment the steam constructor, e.g. `model = MySG(0.0, 2030e3, 1.0, 1.43, 1040.0)`.
3. Set `N` and `M` to match the desired panel (`N=3, M=100` or `N=7, M=50`).
4. Rerun the script to solve and plot with the updated configuration.

## Known gaps

### No automatic publication `savefig`

The appendix scripts produce **interactive Plots** windows. They do not write the publication filenames used in the manuscript (`pressureQuasiWaterN7M50.png`, `MachQuasiSteamN3M100.png`, etc.). Saving Figure 10–11 panels still requires manual `savefig` calls into `figs/`.

### Steam/water toggles

Water vs steam, polynomial degree, and mesh resolution are controlled by manual edits in [`AV_SGQuasiNozzle.jl`](AV_SGQuasiNozzle.jl) rather than command-line arguments.

## Limitations

End-to-end reproduction of Figures 10–11 still requires manual `savefig` calls and steam/water parameter toggles in [`AV_SGQuasiNozzle.jl`](AV_SGQuasiNozzle.jl). The committed `.mat` files support validating the analytical reference data and the exact-solution generator.
