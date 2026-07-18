"""Compare native-angular boundary flux diagnostics.

This example builds a small staged native-angular transport scene with
registered boundary faces. It compares the current no-exchange accounting
baseline, the current sponge source-accounting hook, and the first matched
incident-flux damping candidate on x faces. It also compares an xy-face
no-exchange baseline, an xy-face matched-flux variant, and an x-face
matched-flux variant with a swapped local angular frame so both t and p channel
projections are visible in one run. The example reports contract-driven
incident/outgoing boundary-flux diagnostics and the matched-flux t/p channel
projection. It also reports a static oblique counterpropagating t/p case that
compares summed projected matched flux with direct native-channel matched flux.

The result is an acceptance scene for boundary observability, not a finished
absorber or reflection law.
"""

from __future__ import annotations

import argparse

import numpy as np

import fdtd


def _boundary_face_views(contracts):
    """Return array face slices for boundary contracts with axis/side metadata."""

    axis_index = {"x": 0, "y": 1, "z": 2}
    for contract in contracts:
        axis = contract.get("boundary_axis")
        side = contract.get("boundary_side")
        if axis not in axis_index or side not in ("low", "high"):
            continue
        face = [slice(None), slice(None), slice(None), 0]
        face[axis_index[axis]] = 0 if side == "low" else -1
        yield {
            "boundary_name": contract.get("boundary_name"),
            "face_slice": tuple(face),
        }


def build_boundary_flux_scene(
    boundary_mode: str,
    damping: float = 0.5,
):
    """Build a staged native-angular scene with registered boundary faces."""

    if boundary_mode not in (
        "no_exchange",
        "sponge",
        "matched_flux",
        "matched_flux_swapped_frame",
        "no_exchange_xy",
        "matched_flux_xy",
    ):
        raise ValueError(
            "boundary_mode must be 'no_exchange', 'sponge', 'matched_flux', "
            "'matched_flux_swapped_frame', 'no_exchange_xy', or 'matched_flux_xy'"
        )

    grid = fdtd.AetherGrid(
        shape=(5, 5, 5),
        grid_spacing=8.25e6,
        native_angular_transport=True,
        native_angular_collect_sources=True,
        native_angular_update_clocks=True,
    )
    grid.ell_t[:, :, :, 0] = 0.1
    grid.ell_p[:, :, :, 0] = 0.1
    grid.angular_inertia_t[:, :, :, 0] = 1.0
    grid.angular_inertia_p[:, :, :, 0] = 1.0
    grid.angular_momentum_t[:, :, :, 0] = 1.0
    grid.angular_momentum_p[:, :, :, 0] = 1.0
    if boundary_mode == "matched_flux_swapped_frame":
        grid.configure_native_angular_frame(
            e_t=[0.0, 1.0, 0.0],
            e_p=[1.0, 0.0, 0.0],
        )

    if boundary_mode == "no_exchange":
        grid[0, :, :] = fdtd.AetherAngularNoExchangeBoundary(name="x_low")
        grid[-1, :, :] = fdtd.AetherAngularNoExchangeBoundary(name="x_high")
    elif boundary_mode == "no_exchange_xy":
        grid[0, :, :] = fdtd.AetherAngularNoExchangeBoundary(name="x_low")
        grid[-1, :, :] = fdtd.AetherAngularNoExchangeBoundary(name="x_high")
        grid[:, 0, :] = fdtd.AetherAngularNoExchangeBoundary(name="y_low")
        grid[:, -1, :] = fdtd.AetherAngularNoExchangeBoundary(name="y_high")
    elif boundary_mode == "sponge":
        grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
            damping_t=damping,
            damping_p=damping,
            name="x_low",
        )
        grid[-1, :, :] = fdtd.AetherAngularSpongeBoundary(
            damping_t=damping,
            damping_p=damping,
            name="x_high",
        )
    elif boundary_mode in ("matched_flux", "matched_flux_swapped_frame"):
        grid[0, :, :] = fdtd.AetherAngularMatchedFluxBoundary(
            absorption_t=damping,
            absorption_p=damping,
            name="x_low",
        )
        grid[-1, :, :] = fdtd.AetherAngularMatchedFluxBoundary(
            absorption_t=damping,
            absorption_p=damping,
            name="x_high",
        )
    else:
        grid[0, :, :] = fdtd.AetherAngularMatchedFluxBoundary(
            absorption_t=damping,
            absorption_p=damping,
            name="x_low",
        )
        grid[-1, :, :] = fdtd.AetherAngularMatchedFluxBoundary(
            absorption_t=damping,
            absorption_p=damping,
            name="x_high",
        )
        grid[:, 0, :] = fdtd.AetherAngularMatchedFluxBoundary(
            absorption_t=damping,
            absorption_p=damping,
            name="y_low",
        )
        grid[:, -1, :] = fdtd.AetherAngularMatchedFluxBoundary(
            absorption_t=damping,
            absorption_p=damping,
            name="y_high",
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


def run_boundary_flux_scene(
    boundary_mode: str,
    steps: int = 32,
    damping: float = 0.5,
) -> dict[str, float | int | bool | str]:
    """Advance one boundary-flux scene and return summary diagnostics."""

    grid, previous_t, previous_p = build_boundary_flux_scene(
        boundary_mode=boundary_mode,
        damping=damping,
    )
    contracts = grid.validate_native_angular_boundary_contracts(
        require_physical_laws=boundary_mode
        in ("matched_flux", "matched_flux_swapped_frame", "matched_flux_xy"),
    )
    max_incident = 0.0
    max_outgoing = 0.0
    max_abs_net = 0.0
    max_incident_t = 0.0
    max_incident_p = 0.0
    max_outgoing_t = 0.0
    max_outgoing_p = 0.0
    min_exchange_power = 0.0

    for _ in range(steps):
        current_t = np.asarray(grid.angular_momentum_t).copy()
        current_p = np.asarray(grid.angular_momentum_p).copy()
        grid.update_native_angular_torque(
            previous_momentum_t=previous_t,
            previous_momentum_p=previous_p,
            delta=grid.time_step,
        )
        grid.project_native_angular_torque()

        incident = 0.0
        outgoing = 0.0
        net = 0.0
        incident_t = 0.0
        incident_p = 0.0
        outgoing_t = 0.0
        outgoing_p = 0.0
        for flux in grid.collect_native_angular_boundary_channel_fluxes(
            contracts=contracts,
        ):
            incident += float(flux["incident_total"])
            outgoing += float(flux["outgoing_total"])
            net += float(flux["net_flux"])
            incident_t += float(flux["incident_total_t"])
            incident_p += float(flux["incident_total_p"])
            outgoing_t += float(flux["outgoing_total_t"])
            outgoing_p += float(flux["outgoing_total_p"])
        max_incident = max(max_incident, incident)
        max_outgoing = max(max_outgoing, outgoing)
        max_abs_net = max(max_abs_net, abs(net))
        max_incident_t = max(max_incident_t, incident_t)
        max_incident_p = max(max_incident_p, incident_p)
        max_outgoing_t = max(max_outgoing_t, outgoing_t)
        max_outgoing_p = max(max_outgoing_p, outgoing_p)

        source_t, source_p = grid.collect_native_angular_source_terms()
        _, _, exchange_power = grid.evaluate_native_angular_exchange_power(
            source_t=source_t,
            source_p=source_p,
        )
        min_exchange_power = min(min_exchange_power, float(exchange_power))

        grid.step()
        previous_t = current_t
        previous_p = current_p

    _, _, final_energy = grid.evaluate_native_angular_kinetic_energy()
    momentum_t = np.asarray(grid.angular_momentum_t)
    momentum_p = np.asarray(grid.angular_momentum_p)
    boundary_band_momentum = float(
        sum(
            np.sum(momentum_t[contract["face_slice"]] ** 2)
            + np.sum(momentum_p[contract["face_slice"]] ** 2)
            for contract in _boundary_face_views(contracts)
        )
    )
    return {
        "boundary_mode": boundary_mode,
        "steps": grid.time_steps_passed,
        "boundary_contracts": len(contracts),
        "source_accounting_contracts": bool(
            all(contract["scope"] == "source_accounting" for contract in contracts)
        ),
        "physical_boundary_laws": bool(
            any(contract["physical_boundary_law"] for contract in contracts)
        ),
        "finite_native_momentum": bool(
            np.isfinite(momentum_t).all() and np.isfinite(momentum_p).all()
        ),
        "max_boundary_incident_flux": max_incident,
        "max_boundary_outgoing_flux": max_outgoing,
        "max_abs_boundary_net_flux": max_abs_net,
        "max_boundary_incident_flux_t": max_incident_t,
        "max_boundary_incident_flux_p": max_incident_p,
        "max_boundary_outgoing_flux_t": max_outgoing_t,
        "max_boundary_outgoing_flux_p": max_outgoing_p,
        "min_exchange_power": min_exchange_power,
        "final_native_energy": float(final_energy),
        "boundary_band_momentum": boundary_band_momentum,
    }


def evaluate_counterpropagating_direct_channel_case(
    damping: float = 0.5,
) -> dict[str, float | bool]:
    """Compare weighted and direct channel flux for an oblique t/p cancellation."""

    def _build_grid(boundary_class, name):
        grid = fdtd.AetherGrid(shape=(4, 4, 4), grid_spacing=8.25e6)
        grid.angular_inertia_t[:, :, :, 0] = 1.0
        grid.angular_inertia_p[:, :, :, 0] = 1.0
        grid.angular_momentum_t[:, :, :, 0] = 1.0
        grid.angular_momentum_p[:, :, :, 0] = 1.0
        inv_sqrt_2 = 1.0 / np.sqrt(2.0)
        grid.configure_native_angular_frame(
            e_t=[inv_sqrt_2, inv_sqrt_2, 0.0],
            e_p=[inv_sqrt_2, -inv_sqrt_2, 0.0],
        )
        grid.angular_torque_t[0, 0, 0, 0] = -2.0 * np.sqrt(2.0)
        grid.angular_torque_p[0, 0, 0, 0] = 2.0 * np.sqrt(2.0)
        grid.project_native_angular_torque()
        grid[0, :, :] = boundary_class(
            absorption_t=damping,
            absorption_p=damping,
            name=name,
        )
        return grid

    matched_grid = _build_grid(
        fdtd.AetherAngularMatchedFluxBoundary,
        "x_low_counterpropagating_matched",
    )
    direct_grid = _build_grid(
        fdtd.AetherAngularDirectMatchedFluxBoundary,
        "x_low_counterpropagating_direct",
    )

    matched_contracts = matched_grid.validate_native_angular_boundary_contracts(
        require_physical_laws=True,
    )
    direct_grid.validate_native_angular_boundary_contracts(
        require_physical_laws=True,
    )
    weighted = matched_grid.collect_native_angular_boundary_channel_fluxes(
        contracts=matched_contracts,
    )[0]
    direct = matched_grid.collect_native_angular_boundary_direct_channel_fluxes(
        contracts=matched_contracts,
    )[0]
    matched_t, matched_p = matched_grid.collect_native_angular_source_terms()
    _, _, matched_exchange_power = matched_grid.evaluate_native_angular_exchange_power(
        source_t=matched_t,
        source_p=matched_p,
    )
    direct_t, direct_p = direct_grid.collect_native_angular_source_terms()
    _, _, direct_exchange_power = direct_grid.evaluate_native_angular_exchange_power(
        source_t=direct_t,
        source_p=direct_p,
    )
    matched_source_abs_sum = float(
        np.sum(np.abs(matched_t)) + np.sum(np.abs(matched_p))
    )
    direct_source_abs_sum = float(
        np.sum(np.abs(direct_t)) + np.sum(np.abs(direct_p))
    )

    return {
        "counterpropagating_weighted_incident_total": float(
            weighted["incident_total"]
        ),
        "counterpropagating_weighted_outgoing_total": float(
            weighted["outgoing_total"]
        ),
        "counterpropagating_direct_incident_total": float(
            direct["incident_total"]
        ),
        "counterpropagating_direct_outgoing_total": float(
            direct["outgoing_total"]
        ),
        "counterpropagating_direct_outgoing_t": float(
            direct["outgoing_total_t"]
        ),
        "counterpropagating_direct_incident_p": float(
            direct["incident_total_p"]
        ),
        "counterpropagating_matched_source_abs_sum": matched_source_abs_sum,
        "counterpropagating_matched_exchange_power": float(matched_exchange_power),
        "counterpropagating_direct_matched_source_abs_sum": direct_source_abs_sum,
        "counterpropagating_direct_matched_exchange_power": float(
            direct_exchange_power
        ),
        "counterpropagating_direct_exposes_hidden_flux": bool(
            direct["incident_total"] > weighted["incident_total"]
            and direct["outgoing_total"] > weighted["outgoing_total"]
        ),
        "counterpropagating_matched_no_exchange": bool(
            matched_source_abs_sum == 0.0
            and float(matched_exchange_power) == 0.0
        ),
        "counterpropagating_direct_matched_absorbs_hidden_flux": bool(
            direct_source_abs_sum > 0.0
            and float(direct_exchange_power) < 0.0
        ),
    }


def compare_boundary_flux_scenes(
    steps: int = 32,
    damping: float = 0.5,
) -> dict[str, float | int | bool]:
    """Compare x-face and xy-face boundary-flux scenes."""

    no_exchange = run_boundary_flux_scene(
        boundary_mode="no_exchange",
        steps=steps,
        damping=damping,
    )
    sponge = run_boundary_flux_scene(
        boundary_mode="sponge",
        steps=steps,
        damping=damping,
    )
    matched = run_boundary_flux_scene(
        boundary_mode="matched_flux",
        steps=steps,
        damping=damping,
    )
    matched_swapped = run_boundary_flux_scene(
        boundary_mode="matched_flux_swapped_frame",
        steps=steps,
        damping=damping,
    )
    no_exchange_xy = run_boundary_flux_scene(
        boundary_mode="no_exchange_xy",
        steps=steps,
        damping=damping,
    )
    matched_xy = run_boundary_flux_scene(
        boundary_mode="matched_flux_xy",
        steps=steps,
        damping=damping,
    )
    counterpropagating = evaluate_counterpropagating_direct_channel_case(
        damping=damping,
    )
    return {
        "steps": steps,
        "no_exchange_finite": no_exchange["finite_native_momentum"],
        "sponge_finite": sponge["finite_native_momentum"],
        "matched_flux_finite": matched["finite_native_momentum"],
        "matched_flux_swapped_frame_finite": matched_swapped[
            "finite_native_momentum"
        ],
        "no_exchange_xy_finite": no_exchange_xy["finite_native_momentum"],
        "matched_flux_xy_finite": matched_xy["finite_native_momentum"],
        "no_exchange_boundary_contracts": no_exchange["boundary_contracts"],
        "sponge_boundary_contracts": sponge["boundary_contracts"],
        "matched_flux_boundary_contracts": matched["boundary_contracts"],
        "matched_flux_swapped_frame_boundary_contracts": matched_swapped[
            "boundary_contracts"
        ],
        "no_exchange_xy_boundary_contracts": no_exchange_xy["boundary_contracts"],
        "matched_flux_xy_boundary_contracts": matched_xy["boundary_contracts"],
        "no_exchange_physical_boundary_laws": no_exchange["physical_boundary_laws"],
        "sponge_physical_boundary_laws": sponge["physical_boundary_laws"],
        "matched_flux_physical_boundary_laws": matched["physical_boundary_laws"],
        "matched_flux_swapped_frame_physical_boundary_laws": matched_swapped[
            "physical_boundary_laws"
        ],
        "no_exchange_xy_physical_boundary_laws": no_exchange_xy[
            "physical_boundary_laws"
        ],
        "matched_flux_xy_physical_boundary_laws": matched_xy[
            "physical_boundary_laws"
        ],
        "no_exchange_max_incident_flux": no_exchange["max_boundary_incident_flux"],
        "sponge_max_incident_flux": sponge["max_boundary_incident_flux"],
        "matched_flux_max_incident_flux": matched["max_boundary_incident_flux"],
        "matched_flux_max_incident_flux_t": matched[
            "max_boundary_incident_flux_t"
        ],
        "matched_flux_max_incident_flux_p": matched[
            "max_boundary_incident_flux_p"
        ],
        "matched_flux_swapped_frame_max_incident_flux_t": matched_swapped[
            "max_boundary_incident_flux_t"
        ],
        "matched_flux_swapped_frame_max_incident_flux_p": matched_swapped[
            "max_boundary_incident_flux_p"
        ],
        "matched_flux_xy_max_incident_flux_t": matched_xy[
            "max_boundary_incident_flux_t"
        ],
        "matched_flux_xy_max_incident_flux_p": matched_xy[
            "max_boundary_incident_flux_p"
        ],
        "no_exchange_max_outgoing_flux": no_exchange["max_boundary_outgoing_flux"],
        "sponge_max_outgoing_flux": sponge["max_boundary_outgoing_flux"],
        "matched_flux_max_outgoing_flux": matched["max_boundary_outgoing_flux"],
        "matched_flux_max_outgoing_flux_t": matched[
            "max_boundary_outgoing_flux_t"
        ],
        "matched_flux_max_outgoing_flux_p": matched[
            "max_boundary_outgoing_flux_p"
        ],
        "matched_flux_swapped_frame_max_outgoing_flux_t": matched_swapped[
            "max_boundary_outgoing_flux_t"
        ],
        "matched_flux_swapped_frame_max_outgoing_flux_p": matched_swapped[
            "max_boundary_outgoing_flux_p"
        ],
        "matched_flux_xy_max_outgoing_flux_t": matched_xy[
            "max_boundary_outgoing_flux_t"
        ],
        "matched_flux_xy_max_outgoing_flux_p": matched_xy[
            "max_boundary_outgoing_flux_p"
        ],
        "no_exchange_final_energy": no_exchange["final_native_energy"],
        "sponge_final_energy": sponge["final_native_energy"],
        "matched_flux_final_energy": matched["final_native_energy"],
        "matched_flux_swapped_frame_final_energy": matched_swapped[
            "final_native_energy"
        ],
        "no_exchange_xy_final_energy": no_exchange_xy["final_native_energy"],
        "matched_flux_xy_final_energy": matched_xy["final_native_energy"],
        "sponge_energy_no_greater": bool(
            sponge["final_native_energy"] <= no_exchange["final_native_energy"]
        ),
        "matched_flux_energy_no_greater": bool(
            matched["final_native_energy"] <= no_exchange["final_native_energy"]
        ),
        "matched_flux_swapped_frame_energy_no_greater": bool(
            matched_swapped["final_native_energy"]
            <= no_exchange["final_native_energy"]
        ),
        "matched_flux_xy_energy_no_greater": bool(
            matched_xy["final_native_energy"]
            <= no_exchange_xy["final_native_energy"]
        ),
        "sponge_boundary_momentum_no_greater": bool(
            sponge["boundary_band_momentum"] <= no_exchange["boundary_band_momentum"]
        ),
        "matched_flux_boundary_momentum_no_greater": bool(
            matched["boundary_band_momentum"] <= no_exchange["boundary_band_momentum"]
        ),
        "matched_flux_swapped_frame_boundary_momentum_no_greater": bool(
            matched_swapped["boundary_band_momentum"]
            <= no_exchange["boundary_band_momentum"]
        ),
        "matched_flux_xy_boundary_momentum_no_greater": bool(
            matched_xy["boundary_band_momentum"]
            <= no_exchange_xy["boundary_band_momentum"]
        ),
        "sponge_exchange_nonpositive": bool(sponge["min_exchange_power"] <= 0.0),
        "matched_flux_exchange_nonpositive": bool(
            matched["min_exchange_power"] <= 0.0
        ),
        "matched_flux_swapped_frame_exchange_nonpositive": bool(
            matched_swapped["min_exchange_power"] <= 0.0
        ),
        "matched_flux_xy_exchange_nonpositive": bool(
            matched_xy["min_exchange_power"] <= 0.0
        ),
        **counterpropagating,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare native-angular boundary-flux diagnostics.",
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
        help="number of native-angular transport steps",
    )
    parser.add_argument(
        "--damping",
        type=float,
        default=0.5,
        help="sponge damping coefficient for both native angular channels",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fdtd.set_backend(args.backend)

    summary = compare_boundary_flux_scenes(
        steps=args.steps,
        damping=args.damping,
    )

    print("Native-angular boundary-flux acceptance scene")
    for key, value in summary.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
