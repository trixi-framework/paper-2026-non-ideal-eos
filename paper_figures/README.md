# Main-text figure reproduction

Scripts in this folder reproduce **Figures 1–9** of the compiled paper (`paper/main.pdf` in the parent research repository). Run all commands from this directory.

## Environment setup

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Julia **1.12.6** was used to resolve [`Manifest.toml`](Manifest.toml).

Optional Python environment for jet rendering (from repository root):

```bash
pip install numpy matplotlib pyvista
```

## Figure map

| Paper figure | Script(s) | Primary outputs under `figs/` |
|--------------|-----------|-------------------------------|
| **Figure 1** — EPEC vs APEC (van der Waals density wave) | [`generate_EPEC_vdw_figs.jl`](generate_EPEC_vdw_figs.jl) | `EPEC_comparison_refinement_level_{2,3}_pressure.png` (also total-error panels `EPEC_comparison_refinement_level_*.png` and optional `EPEC_comparison_density_profile.png`) |
| **Figure 2** — APEC vs APEC+LxF (Peng–Robinson wave) | [`generate_APEC_LxF_pr_figs.jl`](generate_APEC_LxF_pr_figs.jl) | `APEC_LxF_comparison_pressure_polydeg_{1,3,7}_refinement_level_{5,4,3}.png` |
| **Figure 3** — APEC vs APEC+entropy correction (smooth/sharp wave) | [`generate_APEC_EC_pr_figs.jl`](generate_APEC_EC_pr_figs.jl) | `blending_coeff_smoothwave_polydeg_{3,7}_refinement_level_{4,3}.png` |
| **Figures 4–5** — Transcritical shock tube | [`generate_pr_shock_figs.jl`](generate_pr_shock_figs.jl) | `transcritical_shock_{rho,p}_polydeg_{3,7}.png`, `transcritical_shock_indicator_flux_{central,terashima_etal,central_terashima_etal}.png` |
| **Figures 6–7** — Transcritical mixing layer | [`trixi_transcritical_mixing.jl`](trixi_transcritical_mixing.jl) | `ec_polydeg3_transcritical_mixing_{rho_rho_c,gamma}.png` and `_long.png` variants (see gaps below) |
| **Figures 8–9** — Transcritical jet | [`trixi_transcritical_jet.jl`](trixi_transcritical_jet.jl), [`plot_vtu_fields.py`](plot_vtu_fields.py) | `transcritical_jet_{rho,pressure,gamma,indicator_shock_capturing}_polydeg{3,7}.png` |

### Shared drivers

| File | Role |
|------|------|
| [`trixi_utils.jl`](trixi_utils.jl) | Custom fluxes, error norms, coordinate helpers |
| [`trixi_transcritical_wave.jl`](trixi_transcritical_wave.jl) | 1D periodic wave driver (van der Waals / Peng–Robinson) |
| [`trixi_transcritical_shock.jl`](trixi_transcritical_shock.jl) | 1D transcritical shock-tube driver |
| [`trixi_transcritical_mixing.jl`](trixi_transcritical_mixing.jl) | 2D mixing-layer simulation and quick plots |
| [`trixi_transcritical_jet.jl`](trixi_transcritical_jet.jl) | 2D jet simulation, HDF5/VTK export, preview plot |

## Run commands

```bash
# Figure 1
julia --project=. generate_EPEC_vdw_figs.jl

# Figure 2
julia --project=. generate_APEC_LxF_pr_figs.jl

# Figure 3 (smooth wave; see manual step for sharp panels)
julia --project=. generate_APEC_EC_pr_figs.jl

# Figures 4–5
julia --project=. generate_pr_shock_figs.jl

# Figures 6–7
julia --project=. trixi_transcritical_mixing.jl

# Figure 8 simulation + VTK export (long runtime)
julia --project=. trixi_transcritical_jet.jl

# Figures 8–9 publication PNGs from committed VTU snapshots
python plot_vtu_fields.py
```

## Precomputed data shipped in this folder

| File | Used for |
|------|----------|
| `solution_000039171.vtu` | Jet Figure 8 (`N=3`) point data |
| `solution_000039171_celldata.vtu` | Jet Figure 8 indicator field |
| `solution_000085663.vtu` | Jet Figure 9 (`N=7`) point data |
| `solution_000085663_celldata.vtu` | Jet Figure 9 indicator field |
| `fast-table-float-0256.csv` | ParaView *Fast* colormap for `plot_vtu_fields.py` |

## Known discrepancies and manual steps

These scripts are copied from the research workspace with **numerical defaults preserved**. The README documents gaps rather than silently changing parameters.

### Figure 1 (`generate_EPEC_vdw_figs.jl`)

- The paper text mentions final time **T = 2** for the van der Waals wave; [`trixi_transcritical_wave.jl`](trixi_transcritical_wave.jl) uses `tspan = (0, 1)` for the van der Waals setup.
- Coppola EPEC flux uses **CFL = 0.001** in the script; the paper discusses stability near **CFL = 0.005**.

### Figure 3 (`generate_APEC_EC_pr_figs.jl`)

- Sharp-wave panels (`blending_coeff_sharpwave_*.png`) require uncommenting `initial_condition_transcritical_wave_sharp` at the top of the generator and rerunning.

### Figures 6–7 (`trixi_transcritical_mixing.jl`)

- Default run produces **2t_c** panels (`ec_polydeg3_transcritical_mixing_*.png`).
- **Figure 7 (4t_c)** requires extending `tspan` to `(0, 0.132)` and uncommenting the `_long.png` `savefig` lines in the driver.
- Degree-7 mixing panels are not produced by the default script configuration.

### Figures 8–9 (jet)

- [`trixi_transcritical_jet.jl`](trixi_transcritical_jet.jl) is hardcoded to **N = 3**, **CFL = 0.4**, and mesh **320×160**. The paper reports **CFL = 0.8** and discusses θ scaling; regenerating Figure 9 requires editing `basis`, `cells_per_dimension`, rerunning the simulation, and selecting the VTK snapshot used in `plot_vtu_fields.py`.
- Publication PNGs for both degrees can be regenerated **without** rerunning Trixi using the committed VTU files and `plot_vtu_fields.py`.

### General

- All `generate_*.jl` scripts call `trixi_include(...)`, provided by **Trixi.jl** for rerunning driver scripts with overridden keyword arguments.
- Intermediate mixing/jet preview PNGs may also appear in this directory; publication outputs belong in `figs/`.

## Limitations

Full end-to-end reproduction of every figure was not re-executed during repository packaging. Expensive 2D simulations (mixing, jet) should be treated as long-running jobs. Use the committed jet VTUs when validating only the Python rendering pipeline.
