"""Adapt the classic quick-start scene to the experimental aether grid.

This example intentionally starts from the source / detector / boundary portion
of ``00-quick-start.ipynb`` so the scene still feels like the original package.
Two limitations are made explicit instead of hidden:

1. ``AetherGrid`` does not yet implement object-aware constitutive updates, so
   the dielectric objects from the notebook are omitted here.
2. Maxwell ``PML`` boundaries are only physically meaningful for the Maxwell
   solver. They can be enabled for the aether path for experimentation, but the
   default aether example leaves them off until substrate-aware absorbing
   boundaries exist.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import fdtd
from fdtd.backend import backend as bd


WAVELENGTH = 1550e-9
SPEED_LIGHT = 299_792_458.0
SOURCE_PERIOD = WAVELENGTH / SPEED_LIGHT


def add_quick_start_boundaries(grid) -> None:
    """Attach the same PML layout used in the original quick-start example."""

    grid[0:10, :, :] = fdtd.PML(name="pml_xlow")
    grid[-10:, :, :] = fdtd.PML(name="pml_xhigh")
    grid[:, 0:10, :] = fdtd.PML(name="pml_ylow")
    grid[:, -10:, :] = fdtd.PML(name="pml_yhigh")


def build_maxwell_quick_start():
    """Build the source/detector/boundary part of ``00-quick-start.ipynb``."""

    grid = fdtd.Grid(shape=(25e-6, 15e-6, 1))
    grid[7.5e-6:8.0e-6, 11.8e-6:13.0e-6, 0] = fdtd.LineSource(
        period=SOURCE_PERIOD,
        phase_shift=np.pi / 2,
        name="source",
    )
    grid[12e-6, :, 0] = fdtd.LineDetector(name="detector")
    add_quick_start_boundaries(grid)
    return grid


def build_aether_quick_start(amplitude: float = 1e-30, use_pml: bool = False):
    """Build the analogous quick-start scene for ``AetherGrid``.

    The geometry and registration style mirror the original quick-start scene,
    but the source is now native to the aether state and the default run does
    not apply Maxwell PMLs.
    """

    grid = fdtd.AetherGrid(shape=(25e-6, 15e-6, 1))
    grid[7.5e-6:8.0e-6, 11.8e-6:13.0e-6, 0] = fdtd.AetherLineSource(
        period=SOURCE_PERIOD,
        phase_shift=np.pi / 2,
        amplitude=amplitude,
        name="source",
    )
    grid[12e-6, :, 0] = fdtd.LineDetector(name="detector")
    if use_pml:
        add_quick_start_boundaries(grid)
    return grid


def run_steps(grid, steps: int) -> None:
    """Advance a grid a fixed number of steps without a progress bar."""

    for _ in range(steps):
        grid.step()


def _detector_peak(detector_history) -> float:
    """Return the peak absolute detector sample over the stored history."""

    if not detector_history:
        return 0.0
    return float(max(np.abs(np.asarray(sample)).max() for sample in detector_history))


def summarize_grid(label: str, grid) -> dict[str, float | int | bool]:
    """Collect a compact numerical summary for terminal output."""

    summary: dict[str, float | int | bool] = {
        "steps": grid.time_steps_passed,
        "finite_E": bool(np.isfinite(np.asarray(grid.E)).all()),
        "finite_H": bool(np.isfinite(np.asarray(grid.H)).all()),
        "max_abs_E": float(np.abs(np.asarray(grid.E)).max()),
        "max_abs_H": float(np.abs(np.asarray(grid.H)).max()),
        "detector_peak_E": _detector_peak(grid.detector.E),
        "detector_peak_H": _detector_peak(grid.detector.H),
    }
    if hasattr(grid, "v"):
        summary["finite_v"] = bool(np.isfinite(np.asarray(grid.v)).all())
        summary["max_abs_v"] = float(np.abs(np.asarray(grid.v)).max())

    print(f"\n{label}")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    return summary


def plot_comparison(maxwell_grid, aether_grid, output: Path, show: bool) -> None:
    """Save a small field snapshot for side-by-side inspection."""

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)

    panels = [
        ("Maxwell Ez", np.asarray(maxwell_grid.E[:, :, 0, 2]), "RdBu"),
        ("Maxwell Hy", np.asarray(maxwell_grid.H[:, :, 0, 1]), "RdBu"),
        ("Aether vz", np.asarray(aether_grid.v[:, :, 0, 2]), "magma"),
        ("Aether Hy", np.asarray(aether_grid.H[:, :, 0, 1]), "RdBu"),
    ]

    for ax, (title, field, cmap) in zip(axes.ravel(), panels):
        m = float(np.abs(field).max()) or 1.0
        im = ax.imshow(field.T, origin="lower", cmap=cmap, vmin=-m, vmax=m)
        ax.set_title(title)
        ax.set_axis_off()
        fig.colorbar(im, ax=ax, shrink=0.8)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=150)
    if show:
        plt.show()
    else:
        plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the quick-start scene on Maxwell Grid and AetherGrid.",
    )
    parser.add_argument(
        "--backend",
        default="numpy",
        help="fdtd backend to use, for example numpy, torch, or torch.cuda",
    )
    parser.add_argument(
        "--maxwell-steps",
        type=int,
        default=40,
        help="number of steps for the Maxwell run",
    )
    parser.add_argument(
        "--aether-steps",
        type=int,
        default=8,
        help="number of steps for the aether run; the current prototype grows rapidly for larger values",
    )
    parser.add_argument(
        "--aether-amplitude",
        type=float,
        default=1e-30,
        help="native source amplitude for the aether run",
    )
    parser.add_argument(
        "--aether-use-pml",
        action="store_true",
        help="apply Maxwell PML boundaries to the aether run for experimentation",
    )
    parser.add_argument(
        "--save-figure",
        type=Path,
        default=Path("fdtd_output/aether_quick_start_comparison.png"),
        help="where to save the comparison figure",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="show the matplotlib window instead of only saving the figure",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fdtd.set_backend(args.backend)

    maxwell_grid = build_maxwell_quick_start()
    aether_grid = build_aether_quick_start(
        amplitude=args.aether_amplitude,
        use_pml=args.aether_use_pml,
    )

    run_steps(maxwell_grid, args.maxwell_steps)
    run_steps(aether_grid, args.aether_steps)

    summarize_grid("Maxwell quick-start subset", maxwell_grid)
    summarize_grid("Aether quick-start adaptation", aether_grid)

    plot_comparison(
        maxwell_grid=maxwell_grid,
        aether_grid=aether_grid,
        output=args.save_figure,
        show=args.show,
    )
    print(f"\nSaved comparison figure to {args.save_figure}")


if __name__ == "__main__":
    main()
