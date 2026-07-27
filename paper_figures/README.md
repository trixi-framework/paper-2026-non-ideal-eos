# Main-text figure reproduction

Scripts in this folder reproduce **Figures 1–9** of the paper. All commands should be run from this directory.

## Environment setup

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Julia 1.12.6 was used to resolve `[Manifest.toml](Manifest.toml)`.

Optional Python environment for jet rendering (from repository root):

```bash
pip install numpy matplotlib pyvista
```

## Figure map


| Paper figure                                                       | Script(s)                                                                                              |
| ------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------ |
| **Figure 1** — EPEC vs APEC (van der Waals density wave)           | `[generate_EPEC_vdw_figs.jl](generate_EPEC_vdw_figs.jl)`                                               |
| **Figure 2** — APEC vs APEC+LxF (Peng–Robinson wave)               | `[generate_APEC_LxF_pr_figs.jl](generate_APEC_LxF_pr_figs.jl)`                                         |
| **Figure 3** — APEC vs APEC+entropy correction (smooth/sharp wave) | `[generate_APEC_EC_pr_figs.jl](generate_APEC_EC_pr_figs.jl)` (loops over smooth and sharp initial conditions) |
| **Figures 4–5** — Transcritical shock tube                         | `[generate_pr_shock_figs.jl](generate_pr_shock_figs.jl)`                                               |
| **Figures 6–7** — Transcritical mixing layer                       | `[trixi_transcritical_mixing.jl](trixi_transcritical_mixing.jl)`                                       |
| **Figures 8–9** — Transcritical jet                                | `[trixi_transcritical_jet.jl](trixi_transcritical_jet.jl)`, `[plot_vtu_fields.py](plot_vtu_fields.py)` |




### Shared drivers


| File                                                             | Role                                                    |
| ---------------------------------------------------------------- | ------------------------------------------------------- |
| `[trixi_utils.jl](trixi_utils.jl)`                               | Custom fluxes, error norms, coordinate helpers          |
| `[trixi_transcritical_wave.jl](trixi_transcritical_wave.jl)`     | 1D periodic wave driver (van der Waals / Peng–Robinson) |
| `[trixi_transcritical_shock.jl](trixi_transcritical_shock.jl)`   | 1D transcritical shock-tube driver                      |
| `[trixi_transcritical_mixing.jl](trixi_transcritical_mixing.jl)` | 2D mixing-layer simulation and quick plots              |
| `[trixi_transcritical_jet.jl](trixi_transcritical_jet.jl)`       | 2D jet simulation, HDF5/VTK export, preview plot        |




## Run commands

```bash
# Figure 1
julia --project=. generate_EPEC_vdw_figs.jl

# Figure 2
julia --project=. generate_APEC_LxF_pr_figs.jl

# Figure 3 (smooth and sharp wave panels)
julia --project=. generate_APEC_EC_pr_figs.jl

# Figures 4–5
julia --project=. generate_pr_shock_figs.jl

# Figures 6–7 (long runtime; optionally start Julia with multiple threads by adding the option `--threads=auto` to the command below)
julia --project=. trixi_transcritical_mixing.jl

# Figure 8 simulation + VTK export (long runtime; optionally start Julia with multiple threads by adding the option `--threads=auto` to the command below)
julia --project=. trixi_transcritical_jet.jl

# Figures 8–9 publication PNGs from committed VTU snapshots
python plot_vtu_fields.py
```



## Precomputed data shipped in this folder


| File                              | Used for                                          |
| --------------------------------- | ------------------------------------------------- |
| `solution_000039171.vtu`          | Jet Figure 8 (`N=3`) point data                   |
| `solution_000039171_celldata.vtu` | Jet Figure 8 indicator field                      |
| `solution_000085663.vtu`          | Jet Figure 9 (`N=7`) point data                   |
| `solution_000085663_celldata.vtu` | Jet Figure 9 indicator field                      |
| `fast-table-float-0256.csv`       | ParaView *Fast* colormap for `plot_vtu_fields.py` |




## Known discrepancies and manual steps

Some paper figures require minor modifications of existing scripts. These are documented below.

### Figures 6–7 (`trixi_transcritical_mixing.jl`)

- By default, this script runs until final time **2t_c** to generate Figures `ec_polydeg3_transcritical_mixing_*.png`. Figure 7 requires running to 4t_c, or extending `tspan` to `(0, 0.132)` and uncommenting the `_long.png` `savefig` lines in the driver.



### Figures 8–9 (jet)

- `[trixi_transcritical_jet.jl](trixi_transcritical_jet.jl)` is hardcoded to N = 3, CFL = 0.8, and mesh 320×160, and must be modified to reproduce the N=7 figure. For convenience, publication PNGs for both degrees can be regenerated without rerunning the elixir using the committed VTU files and `plot_vtu_fields.py`.

