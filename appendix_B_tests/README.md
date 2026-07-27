# Appendix B checks

Scripts in this folder check the Appendix B entropy-Hessian factorizations. Run all commands from this directory.

## Environment setup

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Julia **1.12.6** was used to resolve [`Manifest.toml`](Manifest.toml).

## Run commands

```bash
# Test analytical results numerically
julia --project=. test_appendix_B_eos.jl
```
