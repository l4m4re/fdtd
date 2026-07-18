"""Run a minimal native-angular source / detector scene.

This example exercises the opt-in native-angular transport path without
claiming physical wave propagation. It is intended as a small observable
benchmark: a native-angular point source contributes explicit source terms,
the transport step promotes the native angular momentum candidate, clocks are
resynced from positive inertia, and a native-angular detector records the
result.
"""

from __future__ import annotations

import argparse

import numpy as np

import fdtd


def build_native_angular_observable_scene(
    shape: tuple[int, int, int] = (4, 4, 4),
    grid_spacing: float = 1.0e6,
    amplitude_t: float = 0.5,
    amplitude_p: float = 1.25,
):
    """Build a tiny opt-in native-angular scene."""

    grid = fdtd.AetherGrid(
        shape=shape,
        grid_spacing=grid_spacing,
        native_angular_transport=True,
        native_angular_update_clocks=True,
    )
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0

    source_position = tuple(axis // 2 - 1 for axis in shape)
    grid[source_position] = fdtd.AetherNativeAngularPointSource(
        amplitude_t=amplitude_t,
        amplitude_p=amplitude_p,
        phase_shift=np.pi / 2,
        name="native_source",
    )
    x, y, z = source_position
    grid[x:x, y:y, z:z] = fdtd.AetherNativeAngularDetector(
        name="native_probe",
        fields=(
            "angular_momentum_t",
            "angular_momentum_p",
            "omega_t",
            "omega_p",
            "native_angular_source_t",
            "native_angular_source_p",
            "native_angular_transport_power_t",
            "native_angular_transport_power_p",
        ),
        record_energy=True,
    )
    return grid


def run_scene(
    steps: int = 8,
    amplitude_t: float = 0.5,
    amplitude_p: float = 1.25,
) -> dict[str, float | int | bool]:
    """Advance the observable scene and return a compact numeric summary."""

    grid = build_native_angular_observable_scene(
        amplitude_t=amplitude_t,
        amplitude_p=amplitude_p,
    )

    for _ in range(steps):
        grid.step()

    energy_samples = [float(np.asarray(value)) for value in grid.native_probe.energy]
    summary: dict[str, float | int | bool] = {
        "steps": grid.time_steps_passed,
        "detector_samples": len(grid.native_probe.energy),
        "finite_native_momentum": bool(
            np.isfinite(np.asarray(grid.angular_momentum_t)).all()
            and np.isfinite(np.asarray(grid.angular_momentum_p)).all()
        ),
        "finite_native_clocks": bool(
            np.isfinite(np.asarray(grid.omega_t)).all()
            and np.isfinite(np.asarray(grid.omega_p)).all()
        ),
        "linear_v_quiet": bool(np.allclose(np.asarray(grid.linear_v), 0.0)),
        "max_abs_momentum_t": float(np.abs(np.asarray(grid.angular_momentum_t)).max()),
        "max_abs_momentum_p": float(np.abs(np.asarray(grid.angular_momentum_p)).max()),
        "final_detected_energy": energy_samples[-1] if energy_samples else 0.0,
    }
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a small opt-in native-angular observable scene.",
    )
    parser.add_argument(
        "--backend",
        default="numpy",
        help="fdtd backend to use, for example numpy, torch, or torch.cuda",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=8,
        help="number of native-angular transport steps",
    )
    parser.add_argument(
        "--amplitude-t",
        type=float,
        default=0.5,
        help="toroidal native-angular source amplitude",
    )
    parser.add_argument(
        "--amplitude-p",
        type=float,
        default=1.25,
        help="poloidal native-angular source amplitude",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fdtd.set_backend(args.backend)

    summary = run_scene(
        steps=args.steps,
        amplitude_t=args.amplitude_t,
        amplitude_p=args.amplitude_p,
    )

    print("Native-angular observable scene")
    for key, value in summary.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
