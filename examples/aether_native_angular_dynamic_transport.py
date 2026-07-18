"""Run a staged dynamic native-angular transport scene.

This example is the runnable counterpart of the dynamic Phase 2 transport
tests. It derives native torque from the previous promoted momentum layer,
projects that torque into ``native_angular_tau``, and lets the opt-in
``AetherGrid.step()`` lifecycle promote the next native momentum state.

Optional flags can enable the current sponge source-accounting hook and the
current replacement feedback bridge into ``linear_a``. These remain controlled
experiments on the collocated bridge geometry, not a finished native angular
wave solver.
"""

from __future__ import annotations

import argparse

import numpy as np

import fdtd


def build_dynamic_transport_scene(
    use_sponge: bool = False,
    use_feedback: bool = False,
    feedback_scale: float = 1.0e-26,
):
    """Build the current staged dynamic native-angular transport scene."""

    grid = fdtd.AetherGrid(
        shape=(5, 5, 5),
        grid_spacing=8.25e6,
        native_angular_transport=True,
        native_angular_collect_sources=use_sponge,
        native_angular_update_clocks=True,
        native_angular_feedback_mode="replace" if use_feedback else None,
        native_angular_feedback_scale=feedback_scale,
    )
    grid.ell_t[:, :, :, 0] = 0.1
    grid.ell_p[:, :, :, 0] = 0.1
    grid.angular_inertia_t[:, :, :, 0] = 1.0
    grid.angular_inertia_p[:, :, :, 0] = 1.0
    grid.angular_momentum_t[:, :, :, 0] = 1.0
    grid.angular_momentum_p[:, :, :, 0] = 1.0
    grid[2:2, 2:2, 2:2] = fdtd.AetherNativeAngularDetector(
        name="native_probe",
        fields=(
            "angular_momentum_t",
            "angular_momentum_p",
            "omega_t",
            "omega_p",
            "native_angular_transport_power_t",
            "native_angular_transport_power_p",
        ),
        record_energy=True,
    )

    if use_sponge:
        grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
            damping_t=0.5,
            damping_p=0.5,
        )
        grid[-1, :, :] = fdtd.AetherAngularSpongeBoundary(
            damping_t=0.5,
            damping_p=0.5,
        )

    initial_rate_t = np.zeros((5, 5, 5, 1))
    initial_rate_p = np.zeros((5, 5, 5, 1))
    initial_rate_t[2, 2, 2, 0] = 0.2
    initial_rate_p[2, 2, 2, 0] = -0.1
    previous_t = np.asarray(grid.angular_momentum_t).copy()
    previous_p = np.asarray(grid.angular_momentum_p).copy()
    previous_t -= grid.time_step * initial_rate_t
    previous_p -= grid.time_step * initial_rate_p
    return grid, previous_t, previous_p


def run_dynamic_transport_scene(
    steps: int = 32,
    use_sponge: bool = False,
    use_feedback: bool = False,
    feedback_scale: float = 1.0e-26,
    evaluate_charge: bool = False,
    charge_reduction: str = "additive",
    normalization_area: float = 1.0,
    angular_normalization: float = 1.0,
    loop_measure_t: float = 1.0,
    loop_measure_p: float = 1.0,
) -> dict[str, float | int | bool]:
    """Advance the dynamic scene and return a compact numeric summary."""

    grid, previous_t, previous_p = build_dynamic_transport_scene(
        use_sponge=use_sponge,
        use_feedback=use_feedback,
        feedback_scale=feedback_scale,
    )
    initial_t = np.asarray(grid.angular_momentum_t).copy()
    initial_p = np.asarray(grid.angular_momentum_p).copy()
    max_excursion = 0.0
    max_linear_v = 0.0
    max_response = 0.0
    nonzero_transport_seen = False
    min_exchange_power = 0.0
    max_linear_charge = 0.0
    max_native_charge_t = 0.0
    max_native_charge_p = 0.0
    max_native_charge_reduction = 0.0
    boundary_contracts = grid.validate_native_angular_boundary_contracts()
    max_boundary_incident_flux = 0.0
    max_boundary_outgoing_flux = 0.0
    max_abs_boundary_net_flux = 0.0

    for _ in range(steps):
        current_t = np.asarray(grid.angular_momentum_t).copy()
        current_p = np.asarray(grid.angular_momentum_p).copy()
        grid.update_native_angular_torque(
            previous_momentum_t=previous_t,
            previous_momentum_p=previous_p,
            delta=grid.time_step,
        )
        grid.project_native_angular_torque()
        metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
        nonzero_transport_seen = nonzero_transport_seen or np.any(
            np.asarray(metric_t)
        ) or np.any(np.asarray(metric_p))
        boundary_incident = 0.0
        boundary_outgoing = 0.0
        boundary_net = 0.0
        for flux in grid.collect_native_angular_boundary_fluxes(
            contracts=boundary_contracts
        ):
            boundary_incident += float(flux["incident_total"])
            boundary_outgoing += float(flux["outgoing_total"])
            boundary_net += float(flux["net_flux"])
        max_boundary_incident_flux = max(
            max_boundary_incident_flux,
            boundary_incident,
        )
        max_boundary_outgoing_flux = max(
            max_boundary_outgoing_flux,
            boundary_outgoing,
        )
        max_abs_boundary_net_flux = max(
            max_abs_boundary_net_flux,
            abs(boundary_net),
        )

        if use_sponge:
            source_t, source_p = grid.collect_native_angular_source_terms()
            _, _, exchange_power = grid.evaluate_native_angular_exchange_power(
                source_t=source_t,
                source_p=source_p,
            )
            min_exchange_power = min(min_exchange_power, float(exchange_power))

        grid.step()

        momentum_t = np.asarray(grid.angular_momentum_t)
        momentum_p = np.asarray(grid.angular_momentum_p)
        linear_v = np.asarray(grid.linear_v)
        response = np.asarray(grid.native_angular_linear_response)
        max_excursion = max(
            max_excursion,
            float(np.max(np.abs(momentum_t - initial_t))),
            float(np.max(np.abs(momentum_p - initial_p))),
        )
        max_linear_v = max(max_linear_v, float(np.max(np.abs(linear_v))))
        max_response = max(max_response, float(np.max(np.abs(response))))
        previous_t = current_t
        previous_p = current_p

    if evaluate_charge:
        linear_charge = grid.evaluate_linear_charge_flux_candidate(
            delta_t=grid.time_step,
            area_measure=1.0,
            normalization_area=normalization_area,
            sign=1.0,
        )
        native_charge_t, native_charge_p = (
            grid.evaluate_native_angular_charge_flux_candidates(
                delta_t=grid.time_step,
                loop_measure_t=loop_measure_t,
                loop_measure_p=loop_measure_p,
                angular_normalization=angular_normalization,
                sign_t=1.0,
                sign_p=1.0,
            )
        )
        native_charge_reduction = grid.reduce_native_angular_charge_flux(
            mode=charge_reduction,
        )
        max_linear_charge = float(np.max(np.abs(np.asarray(linear_charge))))
        max_native_charge_t = float(np.max(np.abs(np.asarray(native_charge_t))))
        max_native_charge_p = float(np.max(np.abs(np.asarray(native_charge_p))))
        max_native_charge_reduction = float(
            np.max(np.abs(np.asarray(native_charge_reduction)))
        )

    _, _, final_energy = grid.evaluate_native_angular_kinetic_energy()
    summary: dict[str, float | int | bool] = {
        "steps": grid.time_steps_passed,
        "sponge_enabled": use_sponge,
        "feedback_enabled": use_feedback,
        "detector_samples": len(grid.native_probe.energy),
        "finite_native_momentum": bool(
            np.isfinite(np.asarray(grid.angular_momentum_t)).all()
            and np.isfinite(np.asarray(grid.angular_momentum_p)).all()
        ),
        "finite_linear_v": bool(np.isfinite(np.asarray(grid.linear_v)).all()),
        "nonzero_transport_seen": bool(nonzero_transport_seen),
        "max_abs_momentum_excursion": max_excursion,
        "max_abs_linear_v": max_linear_v,
        "max_abs_linear_response": max_response,
        "min_exchange_power": min_exchange_power,
        "boundary_contracts": len(boundary_contracts),
        "boundary_contracts_source_accounting": bool(
            all(
                contract.get("scope") == "source_accounting"
                for contract in boundary_contracts
            )
        ),
        "physical_boundary_laws": bool(
            any(
                contract.get("physical_boundary_law")
                for contract in boundary_contracts
            )
        ),
        "max_boundary_incident_flux": max_boundary_incident_flux,
        "max_boundary_outgoing_flux": max_boundary_outgoing_flux,
        "max_abs_boundary_net_flux": max_abs_boundary_net_flux,
        "charge_observables_enabled": evaluate_charge,
        "finite_charge_candidates": bool(
            np.isfinite(np.asarray(grid.linear_charge_flux_candidate)).all()
            and np.isfinite(np.asarray(grid.native_angular_charge_flux_t)).all()
            and np.isfinite(np.asarray(grid.native_angular_charge_flux_p)).all()
            and np.isfinite(np.asarray(grid.native_angular_charge_reduction)).all()
        ),
        "max_abs_linear_charge_candidate": max_linear_charge,
        "max_abs_native_charge_t": max_native_charge_t,
        "max_abs_native_charge_p": max_native_charge_p,
        "max_abs_native_charge_reduction": max_native_charge_reduction,
        "final_native_energy": float(final_energy),
    }
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a staged dynamic native-angular transport scene.",
    )
    parser.add_argument(
        "--backend",
        default="numpy",
        help="fdtd backend to use, for example numpy, torch, or torch.cuda",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=32,
        help="number of dynamic native-angular transport steps",
    )
    parser.add_argument(
        "--sponge",
        action="store_true",
        help="enable the current x-face sponge source-accounting hooks",
    )
    parser.add_argument(
        "--feedback",
        action="store_true",
        help="enable opt-in replacement feedback into the linear sector",
    )
    parser.add_argument(
        "--feedback-scale",
        type=float,
        default=1.0e-26,
        help="scale for the opt-in native-angular linear response",
    )
    parser.add_argument(
        "--charge-observables",
        action="store_true",
        help="evaluate explicit local charge-flux candidates after the run",
    )
    parser.add_argument(
        "--charge-reduction",
        choices=("additive", "geometric_mean"),
        default="additive",
        help="native angular charge candidate reduction to report",
    )
    parser.add_argument(
        "--normalization-area",
        type=float,
        default=1.0,
        help="explicit A0 for the local linear charge-flux candidate",
    )
    parser.add_argument(
        "--angular-normalization",
        type=float,
        default=1.0,
        help="explicit Na for native angular charge-flux candidates",
    )
    parser.add_argument(
        "--loop-measure-t",
        type=float,
        default=1.0,
        help="explicit toroidal dtheta/(2*pi)-style loop measure",
    )
    parser.add_argument(
        "--loop-measure-p",
        type=float,
        default=1.0,
        help="explicit poloidal dtheta/(2*pi)-style loop measure",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fdtd.set_backend(args.backend)

    summary = run_dynamic_transport_scene(
        steps=args.steps,
        use_sponge=args.sponge,
        use_feedback=args.feedback,
        feedback_scale=args.feedback_scale,
        evaluate_charge=args.charge_observables,
        charge_reduction=args.charge_reduction,
        normalization_area=args.normalization_area,
        angular_normalization=args.angular_normalization,
        loop_measure_t=args.loop_measure_t,
        loop_measure_p=args.loop_measure_p,
    )

    print("Dynamic native-angular transport scene")
    for key, value in summary.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
