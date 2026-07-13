# Nodal discontinuous Galerkin methods for non-ideal equations of state

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](LICENSE)

Reproducibility repository for the paper

> *Nodal discontinuous Galerkin methods for non-ideal equations of state: pressure equilibrium preservation and entropy correction*

If you find these results useful, please cite the article when it is published. If you use the implementations provided here, please also cite this repository.

## Abstract

This repository contains Julia and Python scripts to reproduce paper results. The workflows are organized into two directories:

| Directory | Contents |
|-----------|----------|
| [`paper_figures/`](paper_figures/) | Main-text figures (Figures 1–9): van der Waals and Peng–Robinson wave/shock tests, mixing layer, and transcritical jet |
| [`appendix_figures/`](appendix_figures/) | Appendix figures (Figures 10–11): stiffened-gas quasi-nozzle benchmarks |

Each directory has its own Julia environment (`Project.toml` and `Manifest.toml`). 

## Requirements

- Julia v1.12.6 was used to generate paper figures and the committed `Manifest.toml` files.
- Python 3 with **NumPy**, **Matplotlib**, and **PyVista** (jet figure post-processing only)

## Quick start

### Main-text figures

```bash
cd paper_figures
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. generate_EPEC_vdw_figs.jl
```

See [`paper_figures/README.md`](paper_figures/README.md) for the full script-to-figure mapping.

### Appendix figures

```bash
cd appendix_figures
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. AV_SGNozzleExact.jl
```

See [`appendix_figures/README.md`](appendix_figures/README.md) for stiffened-gas workflows and known gaps.

### Jet post-processing (Python)

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy matplotlib pyvista
cd paper_figures
python plot_vtu_fields.py
```

Precomputed VTU snapshots for both polynomial degrees are committed beside `plot_vtu_fields.py`.

## Output locations

Generated PNG files are written to `figs/` inside each workflow directory:

- `paper_figures/figs/`
- `appendix_figures/figs/` (when save paths are added manually; see appendix README)

Simulation checkpoints (`out/`, `*.h5`, `*.jld2`) are gitignored and should be regenerated locally.

## Authors

See the paper manuscript for the full author list.

## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
