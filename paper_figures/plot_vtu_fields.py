#!/usr/bin/env python3
"""Render Trixi VTU fields as publication-style PNGs.

Usage:
    pip install numpy matplotlib pyvista
    python plot_vtu_fields.py

VTU files are expected in the same directory as this script.
"""

from __future__ import annotations

from pathlib import Path

from functools import lru_cache

import numpy as np
import pyvista as pv
from matplotlib.colors import ListedColormap

SCRIPT_DIR = Path(__file__).resolve().parent
FIGS_DIR = SCRIPT_DIR / "figs"
FAST_CMAP_CSV = SCRIPT_DIR / "fast-table-float-0256.csv"

SOLUTIONS = [
    {
        "polydeg": 3,
        "point_vtu": SCRIPT_DIR / "solution_000039171.vtu",
        "cell_vtu": SCRIPT_DIR / "solution_000039171_celldata.vtu",
    },
    {
        "polydeg": 7,
        "point_vtu": SCRIPT_DIR / "solution_000085663.vtu",
        "cell_vtu": SCRIPT_DIR / "solution_000085663_celldata.vtu",
    },
]

FIELD_OUTPUTS = [
    ("rho", "transcritical_jet_rho", "point", (16, 660)),
    ("p", "transcritical_jet_pressure", "point", None),
    ("gamma", "transcritical_jet_gamma", "point", (1.4, 50.0)),
    (
        "indicator_shock_capturing",
        "transcritical_jet_indicator_shock_capturing",
        "cell",
        None,
    ),
]

OUTPUTS = [
    (
        solution["point_vtu"] if association == "point" else solution["cell_vtu"],
        field,
        FIGS_DIR / f"{basename}_polydeg{solution['polydeg']}.png",
        association,
        clim,
    )
    for solution in SOLUTIONS
    for field, basename, association, clim in FIELD_OUTPUTS
]

WINDOW_SIZE = (2400, 1300)
SCALAR_BAR_ARGS = {
    "title": "",
    "vertical": False,
    "position_x": 0.15,
    "position_y": 0.005,
    "width": 0.7,
    "height": 0.05,
    "below_label": "",
    "above_label": "",
    "label_font_size": 60,
}

ZOOM_CLOSEST_OFFSET_RATIO = 0.9


def zoom_closest_to_data(plotter: pv.Plotter, mesh: pv.DataSet) -> None:
    """Match ParaView's 'Zoom closest to data' camera framing."""
    plotter.renderer.ResetCameraScreenSpace(mesh.bounds, ZOOM_CLOSEST_OFFSET_RATIO)
    plotter.reset_camera_clipping_range()

@lru_cache(maxsize=1)
def fast_colormap() -> ListedColormap:
    """ParaView 'Fast' preset (Kenneth Moreland)."""
    table = np.genfromtxt(FAST_CMAP_CSV, delimiter=",", skip_header=1)
    return ListedColormap(table[:, 1:4], name="fast")


def fast_lookup_table(vmin: float, vmax: float) -> pv.LookupTable:
    lut = pv.LookupTable()
    lut.cmap = fast_colormap()
    lut.scalar_range = (vmin, vmax)
    return lut


def field_limits(mesh: pv.DataSet, field: str, association: str) -> tuple[float, float]:
    if association == "point":
        data = mesh.point_data[field]
    elif association == "cell":
        data = mesh.cell_data[field]
    else:
        raise ValueError(f"Unknown association: {association!r}")
    return float(np.min(data)), float(np.max(data))


def validate_field(mesh: pv.DataSet, field: str, association: str) -> None:
    container = mesh.point_data if association == "point" else mesh.cell_data
    if field not in container:
        available = ", ".join(container.keys()) or "(none)"
        raise KeyError(
            f"Field {field!r} not found in {association} data. Available: {available}"
        )


def save_field_png(
    vtu_path: Path,
    field: str,
    out_path: Path,
    association: str,
    clim: tuple[float, float] | None = None,
) -> None:
    if not vtu_path.is_file():
        raise FileNotFoundError(f"VTU file not found: {vtu_path}")

    mesh = pv.read(vtu_path)
    validate_field(mesh, field, association)
    if clim is None:
        clim = field_limits(mesh, field, association)
    vmin, vmax = clim
    lut = fast_lookup_table(vmin, vmax)

    out_path.parent.mkdir(parents=True, exist_ok=True)

    plotter = pv.Plotter(off_screen=True, window_size=WINDOW_SIZE)
    plotter.set_background("white")
    plotter.add_mesh(
        mesh,
        scalars=field,
        preference=association,
        lighting=False,
        show_scalar_bar=False,
        cmap=lut,
        clim=(vmin, vmax),
    )
    plotter.view_xy()
    plotter.enable_parallel_projection()
    plotter.show_axes = False
    plotter.remove_all_lights()
    zoom_closest_to_data(plotter, mesh)
    plotter.add_scalar_bar(**SCALAR_BAR_ARGS)
    plotter.screenshot(str(out_path))
    plotter.close()


def main() -> None:
    for vtu_path, field, out_path, association, clim in OUTPUTS:
        print(f"Rendering {field} from {vtu_path.name} -> {out_path.relative_to(SCRIPT_DIR)}")
        save_field_png(vtu_path, field, out_path, association, clim)
    print("Done.")


if __name__ == "__main__":
    main()
