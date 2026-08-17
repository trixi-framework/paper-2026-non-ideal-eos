# Nodal discontinuous Galerkin methods for non-ideal equations of state

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21620203.svg)](https://zenodo.org/doi/10.5281/zenodo.21620203)

Reproducibility repository for the paper

> *Nodal discontinuous Galerkin methods for non-ideal equations of state: pressure equilibrium preservation and entropy correction*

If you find these results useful, please cite the article

```bibtex
@online{chan2026nodal,
  title={Nodal discontinuous {G}alerkin methods for non-ideal equations of state: pressure equilibrium preservation and entropy correction},
  author={Chan, Jesse and Ranocha, Hendrik and Park, Raymond and Lampert, Joshua and Ching, Eric and Edoh, Ayaboe},
  year={2026},
  month={8},
  eprint={2608.14506},
  eprinttype={arxiv},
  eprintclass={math.NA}
}
```

If you use the implementations provided here, please also cite this repository as

```bibtex
@misc{chan2026nodalRepro,
  title={Reproducibility repository for
         "{N}odal discontinuous {G}alerkin methods for non-ideal equations of state: pressure equilibrium preservation and entropy correction"},
  author={Chan, Jesse and Ranocha, Hendrik and Park, Raymond and Lampert, Joshua and Ching, Eric and Edoh, Ayaboe},
  year={2026},
  howpublished={\url{https://github.com/trixi-framework/paper-2026-non-ideal-eos}},
  doi={10.5281/zenodo.21620203}
}
```


## Abstract

Structure-preserving discontinuous Galerkin (DG) methods typically improve the robustness of high order numerical methods for simulations of real fluids.
In addition to conservation, key structures include the preservation of pressure equilibrium and satisfaction of at least one entropy inequality.
In this work, we investigate conservative discretizations using exactly pressure equilibrium conserving (EPEC) and approximately pressure equilibrium conserving (APEC) flux differencing DG formulations, as well as entropy stable formulations through the use of minimally dissipative corrections for non-ideal equations of state (EOS).

We introduce an analysis of EPEC schemes, and compare the performance of two EPEC fluxes.
We also analyze APEC DG schemes and show that the incorporation of dissipative interface penalization terms does not significantly increase pressure equilibrium errors, especially at higher orders of approximation.
Finally, we observe that when combined with APEC flux differencing formulations, entropy correction improves robustness for under-resolved solutions and long-time simulations.


## Numerical experiments

This repository contains Julia and Python scripts to reproduce paper results. The workflows are organized into the following directories:

| Directory | Contents |
|-----------|----------|
| [`paper_figures/`](paper_figures/) | Main-text figures (Figures 1–9): van der Waals and Peng–Robinson wave/shock tests, mixing layer, and transcritical jet |
| [`appendix_figures/`](appendix_figures/) | Appendix figures (Figures 10–11): stiffened-gas quasi-nozzle benchmarks |
| [`appendix_B_tests/`](appendix_B_tests/) | Appendix B verification: analytic entropy-Hessian factorizations for van der Waals and Peng–Robinson (no figures) |

Each directory has its own Julia environment (`Project.toml` and `Manifest.toml`).

## Requirements

- Julia v1.12.6 was used to generate paper figures and the committed `Manifest.toml` files.
- Python 3 with **NumPy**, **Matplotlib**, and **PyVista** (jet figure post-processing only)

## Quick start

### Main-text figures

```bash
cd paper_figures
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. --threads=auto run_generate_figs.jl   # Figures 1–5 (all generate_*.jl)
```

See [`paper_figures/README.md`](paper_figures/README.md) for per-figure commands (including mixing/jet) and the full script-to-figure mapping.

### Appendix figures

```bash
cd appendix_figures
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. --threads=auto run_generate_figs.jl   # exact .mat + default water quasi-nozzle
```

See [`appendix_figures/README.md`](appendix_figures/README.md) for stiffened-gas workflows, other Figure 10–11 panels, and known gaps.

### Appendix B EOS tests

Verifies the Appendix B entropy-Hessian factorizations \(\partial\mathbf{u}/\partial\mathbf{q}\) and \(\partial\mathbf{q}/\partial\mathbf{v}\) for Trixi’s van der Waals and Peng–Robinson EOS in 1D, including `cons2entropy ≈ ∇S`, positive-definiteness of the entropy Hessian, and flux symmetrization. Uses Trixi.jl only (no Clapeyron).

```bash
cd appendix_B_tests
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test_appendix_B_eos.jl
```

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
