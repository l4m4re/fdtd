"""Tests for the experimental aether integration path."""

import os
from pathlib import Path
import subprocess
import sys
from math import pi

import numpy as np

import fdtd
from fdtd.aethergrid import eta
from fdtd.operators import angular_to_linear_bridge, div


class NativeAngularExchangeProbe:
    """Minimal scene element exposing explicit native angular exchange terms."""

    def __init__(self, source_t, source_p, registry):
        self.source_t = source_t
        self.source_p = source_p
        self.registry = registry

    def _register_grid(self, grid, x, y, z):
        self.grid = grid
        self.x = x
        self.y = y
        self.z = z
        getattr(grid, self.registry).append(self)

    def native_angular_source_terms(self):
        return self.source_t, self.source_p


class NativeAngularBoundaryContractProbe:
    """Minimal boundary-like scene element exposing contract metadata."""

    def __init__(self, contract):
        self.contract = contract
        self.name = None

    def _register_grid(self, grid, x, y, z):
        self.grid = grid
        self.x = x
        self.y = y
        self.z = z
        grid.boundaries.append(self)

    def native_angular_boundary_contract(self):
        return self.contract


def test_aethergrid_is_exported():
    grid = fdtd.AetherGrid(shape=(5, 5, 5))
    assert isinstance(grid, fdtd.AetherGrid)


def test_aethergrid_step_advances_time_and_keeps_shapes():
    grid = fdtd.AetherGrid(shape=(5, 5, 5))

    grid.step()

    assert grid.time_steps_passed == 1
    assert grid.v.shape == (5, 5, 5, 3)
    assert grid.p.shape == (5, 5, 5, 1)
    assert grid.yank.shape == (5, 5, 5, 3)
    assert grid.linear_v.shape == (5, 5, 5, 3)
    assert grid.angular_omega.shape == (5, 5, 5, 3)


def test_aethergrid_step_leaves_native_angular_opt_in_path_disabled_by_default():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
    grid[1, 1, 1] = NativeAngularExchangeProbe(
        np.asarray(metric_t).copy(),
        np.asarray(metric_p).copy(),
        "sources",
    )
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()
    initial_linear_a = np.asarray(grid.linear_a).copy()
    grid.native_angular_source_t *= 0.0
    grid.native_angular_source_p *= 0.0

    grid.step()

    np.testing.assert_allclose(np.asarray(grid.angular_momentum_t), initial_momentum_t)
    np.testing.assert_allclose(np.asarray(grid.angular_momentum_p), initial_momentum_p)
    np.testing.assert_allclose(np.asarray(grid.native_angular_source_t), 0.0)
    np.testing.assert_allclose(np.asarray(grid.native_angular_source_p), 0.0)
    np.testing.assert_allclose(np.asarray(grid.linear_a), initial_linear_a)
    assert grid.time_steps_passed == 1


def test_aethergrid_step_runs_opt_in_native_angular_transport():
    grid = fdtd.AetherGrid(
        shape=(4, 4, 4),
        native_angular_transport=True,
    )
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
    metric_t = np.asarray(metric_t).copy()
    metric_p = np.asarray(metric_p).copy()
    grid[1, 1, 1] = NativeAngularExchangeProbe(metric_t, metric_p, "sources")
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()

    grid.step()

    np.testing.assert_allclose(np.asarray(grid.angular_momentum_t), initial_momentum_t)
    np.testing.assert_allclose(np.asarray(grid.angular_momentum_p), initial_momentum_p)
    np.testing.assert_allclose(np.asarray(grid.native_angular_source_t), metric_t)
    np.testing.assert_allclose(np.asarray(grid.native_angular_source_p), metric_p)
    assert grid.time_steps_passed == 1
    assert not np.any(np.asarray(grid.linear_v))


def test_aethergrid_step_native_angular_transport_no_drift_64_step():
    grid = fdtd.AetherGrid(
        shape=(4, 4, 4),
        native_angular_transport=True,
    )
    steps = 64
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
    metric_t = np.asarray(metric_t).copy()
    metric_p = np.asarray(metric_p).copy()
    grid[1, 1, 1] = NativeAngularExchangeProbe(metric_t, metric_p, "sources")
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()

    for _ in range(steps):
        grid.step()
        np.testing.assert_allclose(
            np.asarray(grid.angular_momentum_t),
            initial_momentum_t,
        )
        np.testing.assert_allclose(
            np.asarray(grid.angular_momentum_p),
            initial_momentum_p,
        )
        assert np.all(np.isfinite(np.asarray(grid.angular_momentum_t)))
        assert np.all(np.isfinite(np.asarray(grid.angular_momentum_p)))

    np.testing.assert_allclose(np.asarray(grid.native_angular_source_t), metric_t)
    np.testing.assert_allclose(np.asarray(grid.native_angular_source_p), metric_p)
    assert grid.time_steps_passed == steps
    assert not np.any(np.asarray(grid.linear_v))


def test_aethergrid_step_runs_opt_in_native_angular_feedback():
    grid = fdtd.AetherGrid(
        shape=(4, 4, 4),
        native_angular_feedback_mode="replace",
        native_angular_feedback_scale=0.5,
    )
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    initial_angular_A = np.asarray(grid.angular_A).copy()

    expected_response = 0.5 * angular_to_linear_bridge(
        grid.native_angular_tau,
        metric_length=grid.angular_metric_length,
    )
    initial_native_tau = np.asarray(grid.native_angular_tau).copy()
    grid.step()

    np.testing.assert_allclose(
        np.asarray(grid.native_angular_linear_response),
        np.asarray(expected_response),
    )
    np.testing.assert_allclose(np.asarray(grid.native_angular_tau), initial_native_tau)
    np.testing.assert_allclose(np.asarray(grid.linear_a), np.asarray(expected_response))
    np.testing.assert_allclose(np.asarray(grid.angular_A), initial_angular_A)
    np.testing.assert_allclose(
        np.asarray(grid.linear_v),
        grid.time_step * np.asarray(expected_response)
        + 0.5 * grid.time_step**2 * np.asarray(grid.linear_j),
    )
    assert np.any(np.asarray(grid.linear_v))
    assert grid.time_steps_passed == 1


def test_aethergrid_step_native_angular_feedback_replace_64_step_is_bounded():
    grid = fdtd.AetherGrid(
        shape=(4, 4, 4),
        native_angular_feedback_mode="replace",
        native_angular_feedback_scale=0.5,
    )
    steps = 64
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    expected_response = 0.5 * angular_to_linear_bridge(
        grid.native_angular_tau,
        metric_length=grid.angular_metric_length,
    )
    initial_native_tau = np.asarray(grid.native_angular_tau).copy()

    grid.step()
    first_increment = np.asarray(grid.linear_v).copy()
    np.testing.assert_allclose(
        np.asarray(grid.native_angular_linear_response),
        np.asarray(expected_response),
    )
    assert np.any(first_increment)

    for step in range(1, steps):
        grid.step()
        np.testing.assert_allclose(
            np.asarray(grid.native_angular_linear_response),
            np.asarray(expected_response),
        )
        np.testing.assert_allclose(
            np.asarray(grid.linear_v),
            (step + 1) * first_increment,
        )
        np.testing.assert_allclose(np.asarray(grid.linear_a), np.asarray(expected_response))
        np.testing.assert_allclose(np.asarray(grid.native_angular_tau), initial_native_tau)
        assert np.all(np.isfinite(np.asarray(grid.linear_v)))
        assert np.all(np.isfinite(np.asarray(grid.angular_A)))
        assert np.all(np.isfinite(np.asarray(grid.angular_tau)))

    assert grid.time_steps_passed == steps


def test_aethergrid_rejects_invalid_native_angular_feedback_mode():
    try:
        fdtd.AetherGrid(
            shape=(2, 2, 2),
            native_angular_feedback_mode="invalid",
        )
    except ValueError as exc:
        assert "native_angular_feedback_mode" in str(exc)
    else:
        raise AssertionError("expected native feedback mode validation")


def test_aethergrid_exposes_split_sector_aliases():
    grid = fdtd.AetherGrid(shape=(5, 5, 5))

    assert grid.v is grid.linear_v
    assert grid.p is grid.linear_p
    assert grid.E is grid.linear_E
    assert grid.omega is grid.angular_omega
    assert grid.H is grid.angular_H
    assert grid.A is grid.angular_A
    assert grid.linear_grid_spacing == grid.grid_spacing
    assert grid.angular_grid_spacing == grid.grid_spacing
    assert grid.angular_metric_length.shape == (5, 5, 5, 1)
    assert np.all(np.asarray(grid.angular_metric_length) == grid.angular_grid_spacing)
    assert grid.omega_t.shape == (5, 5, 5, 1)
    assert grid.omega_p.shape == (5, 5, 5, 1)
    assert grid.angular_gamma.shape == (5, 5, 5, 1)
    assert grid.angular_clock_lambda.shape == (5, 5, 5, 1)
    assert grid.angular_e_t.shape == (5, 5, 5, 3)
    assert grid.angular_e_p.shape == (5, 5, 5, 3)
    assert grid.ell_t.shape == (5, 5, 5, 1)
    assert grid.ell_p.shape == (5, 5, 5, 1)
    assert grid.angular_inertia_t.shape == (5, 5, 5, 1)
    assert grid.angular_inertia_p.shape == (5, 5, 5, 1)
    assert grid.angular_momentum_t.shape == (5, 5, 5, 1)
    assert grid.angular_momentum_p.shape == (5, 5, 5, 1)
    assert grid.angular_torque_t.shape == (5, 5, 5, 1)
    assert grid.angular_torque_p.shape == (5, 5, 5, 1)
    assert grid.native_angular_tau.shape == (5, 5, 5, 3)
    assert grid.native_angular_tau_divergence.shape == (5, 5, 5, 1)
    assert grid.native_angular_tau_metric_divergence_t.shape == (5, 5, 5, 1)
    assert grid.native_angular_tau_metric_divergence_p.shape == (5, 5, 5, 1)
    assert grid.native_angular_transport_residual_t.shape == (5, 5, 5, 1)
    assert grid.native_angular_transport_residual_p.shape == (5, 5, 5, 1)
    assert grid.native_angular_source_t.shape == (5, 5, 5, 1)
    assert grid.native_angular_source_p.shape == (5, 5, 5, 1)
    assert grid.native_angular_exchange_power_t.shape == (5, 5, 5, 1)
    assert grid.native_angular_exchange_power_p.shape == (5, 5, 5, 1)
    assert grid.native_angular_momentum_rhs_t.shape == (5, 5, 5, 1)
    assert grid.native_angular_momentum_rhs_p.shape == (5, 5, 5, 1)
    assert grid.native_angular_transport_power_t.shape == (5, 5, 5, 1)
    assert grid.native_angular_transport_power_p.shape == (5, 5, 5, 1)
    assert grid.native_angular_momentum_candidate_t.shape == (5, 5, 5, 1)
    assert grid.native_angular_momentum_candidate_p.shape == (5, 5, 5, 1)
    assert grid.native_angular_linear_response.shape == (5, 5, 5, 3)
    assert grid.linear_charge_flux_candidate.shape == (5, 5, 5, 3)
    assert grid.native_angular_charge_flux_t.shape == (5, 5, 5, 1)
    assert grid.native_angular_charge_flux_p.shape == (5, 5, 5, 1)
    assert grid.native_angular_charge_reduction.shape == (5, 5, 5, 1)
    assert grid.native_angular_boundary_direct_normal_flux_t.shape == (5, 5, 5, 1)
    assert grid.native_angular_boundary_direct_normal_flux_p.shape == (5, 5, 5, 1)
    assert grid.native_angular_boundary_direct_incident_flux_t.shape == (
        5,
        5,
        5,
        1,
    )
    assert grid.native_angular_boundary_direct_incident_flux_p.shape == (
        5,
        5,
        5,
        1,
    )
    assert grid.native_angular_boundary_direct_outgoing_flux_t.shape == (
        5,
        5,
        5,
        1,
    )
    assert grid.native_angular_boundary_direct_outgoing_flux_p.shape == (
        5,
        5,
        5,
        1,
    )


def test_aethergrid_exposes_default_native_angular_geometry():
    grid = fdtd.AetherGrid(shape=(3, 3, 3))

    np.testing.assert_allclose(np.asarray(grid.angular_e_t[..., 0]), 1.0)
    np.testing.assert_allclose(np.asarray(grid.angular_e_t[..., 1:]), 0.0)
    np.testing.assert_allclose(np.asarray(grid.angular_e_p[..., 1]), 1.0)
    np.testing.assert_allclose(np.asarray(grid.angular_e_p[..., (0, 2)]), 0.0)
    np.testing.assert_allclose(np.asarray(grid.ell_t), grid.angular_grid_spacing)
    np.testing.assert_allclose(np.asarray(grid.ell_p), grid.angular_grid_spacing)
    assert not np.any(np.asarray(grid.angular_inertia_t))
    assert not np.any(np.asarray(grid.angular_inertia_p))
    assert not np.any(np.asarray(grid.angular_momentum_t))
    assert not np.any(np.asarray(grid.angular_momentum_p))
    assert not np.any(np.asarray(grid.angular_torque_t))
    assert not np.any(np.asarray(grid.angular_torque_p))
    assert not np.any(np.asarray(grid.native_angular_tau))
    assert not np.any(np.asarray(grid.native_angular_tau_divergence))
    assert not np.any(np.asarray(grid.native_angular_tau_metric_divergence_t))
    assert not np.any(np.asarray(grid.native_angular_tau_metric_divergence_p))
    assert not np.any(np.asarray(grid.native_angular_transport_residual_t))
    assert not np.any(np.asarray(grid.native_angular_transport_residual_p))
    assert not np.any(np.asarray(grid.native_angular_source_t))
    assert not np.any(np.asarray(grid.native_angular_source_p))
    assert not np.any(np.asarray(grid.native_angular_exchange_power_t))
    assert not np.any(np.asarray(grid.native_angular_exchange_power_p))
    assert not np.any(np.asarray(grid.native_angular_momentum_rhs_t))
    assert not np.any(np.asarray(grid.native_angular_momentum_rhs_p))
    assert not np.any(np.asarray(grid.native_angular_transport_power_t))
    assert not np.any(np.asarray(grid.native_angular_transport_power_p))
    assert not np.any(np.asarray(grid.native_angular_momentum_candidate_t))
    assert not np.any(np.asarray(grid.native_angular_momentum_candidate_p))
    assert not np.any(np.asarray(grid.native_angular_linear_response))
    assert not np.any(np.asarray(grid.linear_charge_flux_candidate))
    assert not np.any(np.asarray(grid.native_angular_charge_flux_t))
    assert not np.any(np.asarray(grid.native_angular_charge_flux_p))
    assert not np.any(np.asarray(grid.native_angular_charge_reduction))


def test_aethergrid_evaluates_native_angular_clock_benchmark():
    grid = fdtd.AetherGrid(shape=(3, 3, 3))
    delta = 0.2
    grid.omega_t[1, 1, 1, 0] = 1.2
    grid.omega_p[1, 1, 1, 0] = 0.8
    grid.angular_gamma[1, 1, 1, 0] = 0.4

    result = grid.evaluate_angular_clock_benchmark(delta=delta)

    expected = (
        -4.0 / delta**2 * np.sin(delta * 1.2 / 2.0) ** 2
        -4.0 / delta**2 * np.sin(delta * 0.8 / 2.0) ** 2
        +4.0 / delta**2 * np.sinh(delta * 0.4 / 2.0) ** 2
    )
    assert float(result[1, 1, 1, 0]) == np.asarray(
        grid.angular_clock_lambda
    )[1, 1, 1, 0]
    np.testing.assert_allclose(float(result[1, 1, 1, 0]), expected)
    np.testing.assert_allclose(float(grid.theta_t[1, 1, 1, 0]), 1.2 * delta)
    np.testing.assert_allclose(float(grid.theta_p[1, 1, 1, 0]), 0.8 * delta)
    np.testing.assert_allclose(float(grid.angular_chi[1, 1, 1, 0]), 0.4 * delta)


def test_aethergrid_updates_native_angular_inertia_from_metric_lengths():
    grid = fdtd.AetherGrid(shape=(3, 3, 3))
    density = 2.0
    grid.ell_t[1, 1, 1, 0] = 3.0
    grid.ell_p[1, 1, 1, 0] = 4.0

    inertia_t, inertia_p = grid.update_native_angular_inertia(density=density)

    assert float(inertia_t[1, 1, 1, 0]) == 18.0
    assert float(inertia_p[1, 1, 1, 0]) == 32.0
    assert float(grid.angular_inertia_t[1, 1, 1, 0]) == 18.0
    assert float(grid.angular_inertia_p[1, 1, 1, 0]) == 32.0


def test_aethergrid_native_angular_inertia_accepts_density_field():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))
    density = np.ones((2, 2, 2, 1)) * 3.0
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 5.0

    inertia_t, inertia_p = grid.update_native_angular_inertia(density=density)

    np.testing.assert_allclose(np.asarray(inertia_t), 12.0)
    np.testing.assert_allclose(np.asarray(inertia_p), 75.0)


def test_aethergrid_updates_native_angular_momentum_from_inertia_and_clocks():
    grid = fdtd.AetherGrid(shape=(3, 3, 3))
    grid.angular_inertia_t[1, 1, 1, 0] = 6.0
    grid.angular_inertia_p[1, 1, 1, 0] = 8.0
    grid.omega_t[1, 1, 1, 0] = 1.5
    grid.omega_p[1, 1, 1, 0] = 2.5

    momentum_t, momentum_p = grid.update_native_angular_momentum()

    assert float(momentum_t[1, 1, 1, 0]) == 9.0
    assert float(momentum_p[1, 1, 1, 0]) == 20.0
    assert float(grid.angular_momentum_t[1, 1, 1, 0]) == 9.0
    assert float(grid.angular_momentum_p[1, 1, 1, 0]) == 20.0


def test_aethergrid_evaluates_native_angular_kinetic_energy_diagnostic():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 6.0
    grid.angular_momentum_p[:, :, :, 0] = 8.0

    energy_t, energy_p, total_energy = (
        grid.evaluate_native_angular_kinetic_energy()
    )

    np.testing.assert_allclose(np.asarray(energy_t), 9.0)
    np.testing.assert_allclose(np.asarray(energy_p), 8.0)
    assert float(total_energy) == 8 * (9.0 + 8.0)


def test_aethergrid_native_angular_kinetic_energy_requires_positive_inertia():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))
    grid.angular_inertia_t[:, :, :, 0] = 1.0
    grid.angular_inertia_p[:, :, :, 0] = 0.0

    try:
        grid.evaluate_native_angular_kinetic_energy()
    except ValueError as exc:
        assert "inertia must be positive" in str(exc)
    else:
        raise AssertionError("expected positive-inertia validation")


def test_aethergrid_evaluates_native_angular_exchange_power_diagnostic():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))
    source_t = np.zeros((2, 2, 2, 1))
    source_p = np.zeros((2, 2, 2, 1))
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[1, 1, 1, 0] = 6.0
    grid.angular_momentum_p[1, 1, 1, 0] = 8.0
    source_t[1, 1, 1, 0] = 0.5
    source_p[1, 1, 1, 0] = -1.5

    power_t, power_p, total_power = grid.evaluate_native_angular_exchange_power(
        source_t=source_t,
        source_p=source_p,
    )

    np.testing.assert_allclose(float(power_t[1, 1, 1, 0]), 1.5)
    np.testing.assert_allclose(float(power_p[1, 1, 1, 0]), -3.0)
    np.testing.assert_allclose(float(total_power), -1.5)
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_native_angular_exchange_power_requires_positive_inertia():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))
    grid.angular_inertia_t[:, :, :, 0] = 1.0
    grid.angular_inertia_p[:, :, :, 0] = 0.0

    try:
        grid.evaluate_native_angular_exchange_power()
    except ValueError as exc:
        assert "inertia must be positive" in str(exc)
    else:
        raise AssertionError("expected positive-inertia validation")


def test_aethergrid_evaluates_native_angular_transport_power_diagnostic():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    source_t = np.ones((4, 4, 4, 1)) * 0.1
    source_p = np.ones((4, 4, 4, 1)) * 0.2

    power_t, power_p, total_power = grid.evaluate_native_angular_transport_power(
        source_t=source_t,
        source_p=source_p,
    )

    expected_divergence = np.asarray(div(grid.native_angular_tau))
    expected_rhs_t = source_t - 2.0 * expected_divergence
    expected_rhs_p = source_p - 3.0 * expected_divergence
    expected_power_t = expected_rhs_t * 0.75 / 2.0
    expected_power_p = expected_rhs_p * 1.25 / 4.0
    np.testing.assert_allclose(np.asarray(power_t), expected_power_t)
    np.testing.assert_allclose(np.asarray(power_p), expected_power_p)
    np.testing.assert_allclose(
        float(total_power),
        np.sum(expected_power_t + expected_power_p),
    )
    assert power_t is grid.native_angular_transport_power_t
    assert power_p is grid.native_angular_transport_power_p
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_native_angular_transport_power_requires_positive_inertia():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))
    grid.angular_inertia_t[:, :, :, 0] = 1.0
    grid.angular_inertia_p[:, :, :, 0] = 0.0

    try:
        grid.evaluate_native_angular_transport_power()
    except ValueError as exc:
        assert "inertia must be positive" in str(exc)
    else:
        raise AssertionError("expected positive-inertia validation")


def test_aethergrid_updates_native_angular_torque_from_momentum_change():
    grid = fdtd.AetherGrid(shape=(3, 3, 3))
    previous_t = np.zeros((3, 3, 3, 1))
    previous_p = np.zeros((3, 3, 3, 1))
    previous_t[1, 1, 1, 0] = 2.0
    previous_p[1, 1, 1, 0] = 4.0
    grid.angular_momentum_t[1, 1, 1, 0] = 8.0
    grid.angular_momentum_p[1, 1, 1, 0] = 10.0

    torque_t, torque_p = grid.update_native_angular_torque(
        previous_momentum_t=previous_t,
        previous_momentum_p=previous_p,
        delta=0.5,
    )

    assert float(torque_t[1, 1, 1, 0]) == 12.0
    assert float(torque_p[1, 1, 1, 0]) == 12.0
    assert float(grid.angular_torque_t[1, 1, 1, 0]) == 12.0
    assert float(grid.angular_torque_p[1, 1, 1, 0]) == 12.0


def test_aethergrid_native_angular_torque_defaults_previous_momentum_to_zero():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))
    grid.angular_momentum_t[:, :, :, 0] = 3.0
    grid.angular_momentum_p[:, :, :, 0] = 6.0

    torque_t, torque_p = grid.update_native_angular_torque(delta=3.0)

    np.testing.assert_allclose(np.asarray(torque_t), 1.0)
    np.testing.assert_allclose(np.asarray(torque_p), 2.0)


def test_aethergrid_native_angular_torque_rejects_zero_step():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))

    try:
        grid.update_native_angular_torque(delta=0.0)
    except ValueError as exc:
        assert "delta" in str(exc)
    else:
        raise AssertionError("expected ValueError for zero delta")


def test_aethergrid_projects_native_angular_torque_to_local_frame():
    grid = fdtd.AetherGrid(shape=(3, 3, 3))
    grid.angular_torque_t[1, 1, 1, 0] = 2.0
    grid.angular_torque_p[1, 1, 1, 0] = 3.0

    projected = grid.project_native_angular_torque()

    np.testing.assert_allclose(np.asarray(projected[1, 1, 1]), [2.0, 3.0, 0.0])
    assert projected is grid.native_angular_tau
    assert not np.any(np.asarray(grid.angular_tau))


def test_aethergrid_projects_native_angular_torque_with_custom_frame():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))
    grid.angular_e_t[0, 0, 0] = [0.0, 0.0, 1.0]
    grid.angular_e_p[0, 0, 0] = [1.0, 0.0, 0.0]
    grid.angular_torque_t[0, 0, 0, 0] = 2.0
    grid.angular_torque_p[0, 0, 0, 0] = 3.0

    projected = grid.project_native_angular_torque()

    np.testing.assert_allclose(np.asarray(projected[0, 0, 0]), [3.0, 0.0, 2.0])


def test_aethergrid_configures_uniform_native_angular_frame():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))

    e_t, e_p = grid.configure_native_angular_frame(
        e_t=[0.0, 1.0, 0.0],
        e_p=[1.0, 0.0, 0.0],
    )

    np.testing.assert_allclose(np.asarray(e_t[..., 1]), 1.0)
    np.testing.assert_allclose(np.asarray(e_t[..., (0, 2)]), 0.0)
    np.testing.assert_allclose(np.asarray(e_p[..., 0]), 1.0)
    np.testing.assert_allclose(np.asarray(e_p[..., 1:]), 0.0)


def test_aethergrid_configures_field_native_angular_frame():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))
    e_t = np.zeros((2, 2, 2, 3))
    e_p = np.zeros((2, 2, 2, 3))
    e_t[..., 2] = 1.0
    e_p[..., 0] = 1.0

    grid.configure_native_angular_frame(e_t=e_t, e_p=e_p)
    grid.angular_torque_t[0, 0, 0, 0] = 2.0
    grid.angular_torque_p[0, 0, 0, 0] = 3.0

    projected = grid.project_native_angular_torque()

    np.testing.assert_allclose(np.asarray(projected[0, 0, 0]), [3.0, 0.0, 2.0])


def test_aethergrid_native_angular_frame_configuration_validates_inputs():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))

    cases = [
        (
            [1.0, 0.0],
            [0.0, 1.0, 0.0],
            "e_t must have shape",
        ),
        (
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            "nonzero",
        ),
        (
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            "orthogonal",
        ),
    ]
    for e_t, e_p, message in cases:
        try:
            grid.configure_native_angular_frame(e_t=e_t, e_p=e_p)
        except ValueError as exc:
            assert message in str(exc)
        else:
            raise AssertionError("expected frame validation")


def test_aethergrid_evaluates_native_angular_torque_divergence():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0

    diagnostic = grid.evaluate_native_angular_torque_divergence()

    np.testing.assert_allclose(
        np.asarray(diagnostic),
        np.asarray(div(grid.native_angular_tau)),
    )
    assert diagnostic is grid.native_angular_tau_divergence
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_evaluates_metric_weighted_torque_divergence():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0

    metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
    expected_divergence = np.asarray(div(grid.native_angular_tau))

    np.testing.assert_allclose(np.asarray(metric_t), 2.0 * expected_divergence)
    np.testing.assert_allclose(np.asarray(metric_p), 3.0 * expected_divergence)
    assert metric_t is grid.native_angular_tau_metric_divergence_t
    assert metric_p is grid.native_angular_tau_metric_divergence_p
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_uniform_native_angular_torque_has_zero_divergence():
    grid = fdtd.AetherGrid(shape=(5, 5, 5))
    grid.native_angular_tau[:, :, :, 0] = 2.0
    grid.native_angular_tau[:, :, :, 1] = -3.0
    grid.native_angular_tau[:, :, :, 2] = 5.0

    diagnostic = grid.evaluate_native_angular_torque_divergence()
    metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()

    np.testing.assert_allclose(np.asarray(diagnostic), 0.0)
    np.testing.assert_allclose(np.asarray(metric_t), 0.0)
    np.testing.assert_allclose(np.asarray(metric_p), 0.0)
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_linear_native_torque_transport_uses_current_stencil():
    grid = fdtd.AetherGrid(shape=(5, 5, 5))
    x_pattern = np.array([1.0, 2.0, 2.0, 1.0, 0.0])
    y_pattern = np.array([2.0, 4.0, 4.0, 2.0, 0.0])
    z_pattern = np.array([3.0, 6.0, 6.0, 3.0, 0.0])
    expected_divergence = (
        x_pattern[:, None, None]
        + y_pattern[None, :, None]
        + z_pattern[None, None, :]
    )[:, :, :, None]
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_torque_t[:, :, :, 0] = 0.25
    grid.angular_torque_p[:, :, :, 0] = 0.5
    for x in range(5):
        grid.native_angular_tau[x, :, :, 0] = x
    for y in range(5):
        grid.native_angular_tau[:, y, :, 1] = 2.0 * y
    for z in range(5):
        grid.native_angular_tau[:, :, z, 2] = 3.0 * z

    diagnostic = grid.evaluate_native_angular_torque_divergence()
    metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
    source_t = np.asarray(grid.angular_torque_t) + np.asarray(metric_t)
    source_p = np.asarray(grid.angular_torque_p) + np.asarray(metric_p)
    residual_t, residual_p = grid.evaluate_native_angular_transport_residual(
        source_t=source_t,
        source_p=source_p,
    )

    np.testing.assert_allclose(np.asarray(diagnostic), expected_divergence)
    np.testing.assert_allclose(np.asarray(metric_t), 2.0 * expected_divergence)
    np.testing.assert_allclose(np.asarray(metric_p), 3.0 * expected_divergence)
    np.testing.assert_allclose(np.asarray(residual_t), 0.0)
    np.testing.assert_allclose(np.asarray(residual_p), 0.0)
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_candidate_only_native_torque_transport_loop_is_bounded():
    grid = fdtd.AetherGrid(shape=(5, 5, 5))
    delta = 0.01
    steps = 32
    grid.ell_t[:, :, :, 0] = 0.1
    grid.ell_p[:, :, :, 0] = 0.1
    grid.angular_momentum_t[:, :, :, 0] = 1.0
    grid.angular_momentum_p[:, :, :, 0] = 1.0
    initial_t = np.asarray(grid.angular_momentum_t).copy()
    initial_p = np.asarray(grid.angular_momentum_p).copy()
    initial_rate_t = np.zeros((5, 5, 5, 1))
    initial_rate_p = np.zeros((5, 5, 5, 1))
    initial_rate_t[2, 2, 2, 0] = 0.2
    initial_rate_p[2, 2, 2, 0] = -0.1
    previous_t = initial_t - delta * initial_rate_t
    previous_p = initial_p - delta * initial_rate_p
    max_excursion = 0.0
    nonzero_transport_seen = False

    for _ in range(steps):
        grid.update_native_angular_torque(
            previous_momentum_t=previous_t,
            previous_momentum_p=previous_p,
            delta=delta,
        )
        grid.project_native_angular_torque()
        metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
        nonzero_transport_seen = nonzero_transport_seen or np.any(
            np.asarray(metric_t)
        ) or np.any(np.asarray(metric_p))
        candidate_t, candidate_p = grid.predict_native_angular_momentum_step(
            delta=delta,
        )
        candidate_t = np.asarray(candidate_t).copy()
        candidate_p = np.asarray(candidate_p).copy()

        assert np.all(np.isfinite(candidate_t))
        assert np.all(np.isfinite(candidate_p))
        max_excursion = max(
            max_excursion,
            float(np.max(np.abs(candidate_t - initial_t))),
            float(np.max(np.abs(candidate_p - initial_p))),
        )

        # Candidate-only propagated-transport benchmark; production dynamics
        # still do not call or apply the predictor.
        previous_t = np.asarray(grid.angular_momentum_t).copy()
        previous_p = np.asarray(grid.angular_momentum_p).copy()
        grid.angular_momentum_t = grid.native_angular_momentum_candidate_t
        grid.angular_momentum_p = grid.native_angular_momentum_candidate_p

    assert nonzero_transport_seen
    assert 0.0 < max_excursion < 1.0e-3
    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))
    assert not np.any(np.asarray(grid.native_angular_source_t))
    assert not np.any(np.asarray(grid.native_angular_source_p))


def test_aethergrid_step_dynamic_native_torque_transport_loop_is_bounded():
    grid = fdtd.AetherGrid(
        shape=(5, 5, 5),
        grid_spacing=8.25e6,
        native_angular_transport=True,
        native_angular_collect_sources=False,
        native_angular_update_clocks=True,
    )
    steps = 32
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
        ),
        record_energy=True,
    )
    initial_t = np.asarray(grid.angular_momentum_t).copy()
    initial_p = np.asarray(grid.angular_momentum_p).copy()
    initial_rate_t = np.zeros((5, 5, 5, 1))
    initial_rate_p = np.zeros((5, 5, 5, 1))
    initial_rate_t[2, 2, 2, 0] = 0.2
    initial_rate_p[2, 2, 2, 0] = -0.1
    previous_t = initial_t - grid.time_step * initial_rate_t
    previous_p = initial_p - grid.time_step * initial_rate_p
    max_excursion = 0.0
    nonzero_transport_seen = False

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

        grid.step()

        momentum_t = np.asarray(grid.angular_momentum_t)
        momentum_p = np.asarray(grid.angular_momentum_p)
        assert np.all(np.isfinite(momentum_t))
        assert np.all(np.isfinite(momentum_p))
        max_excursion = max(
            max_excursion,
            float(np.max(np.abs(momentum_t - initial_t))),
            float(np.max(np.abs(momentum_p - initial_p))),
        )
        np.testing.assert_allclose(np.asarray(grid.omega_t), momentum_t)
        np.testing.assert_allclose(np.asarray(grid.omega_p), momentum_p)

        previous_t = current_t
        previous_p = current_p

    assert nonzero_transport_seen
    assert 0.0 < max_excursion < 1.0e-3
    assert grid.time_steps_passed == steps
    assert len(grid.native_probe.energy) == steps
    assert len(grid.native_probe.readings["angular_momentum_t"]) == steps
    assert not np.any(np.asarray(grid.native_angular_source_t))
    assert not np.any(np.asarray(grid.native_angular_source_p))
    assert not np.any(np.asarray(grid.linear_v))


def test_aethergrid_evaluates_passive_native_angular_linear_response():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    metric_length = np.ones((4, 4, 4, 1)) * 0.25
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    grid.native_angular_tau[1, 2, 1, 2] = 0.5
    scale = 0.75

    response = grid.evaluate_native_angular_linear_response(
        metric_length=metric_length,
        scale=scale,
    )

    expected = scale * angular_to_linear_bridge(
        grid.native_angular_tau,
        metric_length=metric_length,
    )
    np.testing.assert_allclose(np.asarray(response), np.asarray(expected))
    assert response is grid.native_angular_linear_response
    assert np.any(np.asarray(response))
    assert not np.any(np.asarray(grid.angular_A))
    assert not np.any(np.asarray(grid.linear_a))
    assert grid.time_steps_passed == 0


def test_aethergrid_native_angular_linear_response_is_coupling_diagnostic_only():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.linear_a[:, :, :, 0] = 7.0
    grid.angular_A[:, :, :, 1] = 11.0
    initial_linear_a = np.asarray(grid.linear_a).copy()
    initial_angular_A = np.asarray(grid.angular_A).copy()
    grid.native_angular_tau[1, 1, 1, 0] = 2.0

    response = grid.evaluate_native_angular_linear_response(scale=0.5)

    assert np.any(np.asarray(response))
    np.testing.assert_allclose(np.asarray(grid.linear_a), initial_linear_a)
    np.testing.assert_allclose(np.asarray(grid.angular_A), initial_angular_A)
    assert not np.any(np.asarray(grid.angular_tau))
    assert grid.time_steps_passed == 0


def test_aethergrid_applies_native_angular_linear_response_additively():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    response = np.ones((4, 4, 4, 3)) * 0.25
    grid.linear_a[:, :, :, 0] = 7.0
    grid.angular_A[:, :, :, 1] = 11.0
    initial_angular_A = np.asarray(grid.angular_A).copy()

    updated = grid.apply_native_angular_linear_response(
        response=response,
        mode="add",
    )

    expected = np.zeros((4, 4, 4, 3))
    expected[:, :, :, 0] = 7.0
    expected += response
    np.testing.assert_allclose(np.asarray(updated), expected)
    np.testing.assert_allclose(np.asarray(grid.linear_a), expected)
    np.testing.assert_allclose(np.asarray(grid.native_angular_linear_response), response)
    np.testing.assert_allclose(np.asarray(grid.angular_A), initial_angular_A)
    assert grid.a is grid.linear_a
    assert not np.any(np.asarray(grid.angular_tau))
    assert grid.time_steps_passed == 0


def test_aethergrid_applies_native_angular_linear_response_by_replacement():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    metric_length = np.ones((4, 4, 4, 1)) * 0.25
    grid.linear_a[:, :, :, 0] = 7.0
    grid.angular_A[:, :, :, 1] = 11.0
    initial_angular_A = np.asarray(grid.angular_A).copy()
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    scale = 0.5

    updated = grid.apply_native_angular_linear_response(
        metric_length=metric_length,
        scale=scale,
        mode="replace",
    )

    expected = scale * angular_to_linear_bridge(
        grid.native_angular_tau,
        metric_length=metric_length,
    )
    np.testing.assert_allclose(np.asarray(updated), np.asarray(expected))
    assert updated is grid.linear_a
    assert updated is grid.native_angular_linear_response
    np.testing.assert_allclose(np.asarray(grid.angular_A), initial_angular_A)
    assert not np.any(np.asarray(grid.angular_tau))
    assert grid.time_steps_passed == 0


def test_aethergrid_native_angular_linear_response_rejects_invalid_apply_inputs():
    grid = fdtd.AetherGrid(shape=(3, 3, 3))
    grid.linear_a[:, :, :, 0] = 7.0
    grid.native_angular_linear_response[:, :, :, 1] = 5.0
    initial_linear_a = np.asarray(grid.linear_a).copy()
    initial_response = np.asarray(grid.native_angular_linear_response).copy()

    try:
        grid.apply_native_angular_linear_response(
            response=np.ones((3, 3, 3, 1)),
        )
    except ValueError as exc:
        assert "shape" in str(exc)
    else:
        raise AssertionError("expected response-shape validation")
    np.testing.assert_allclose(np.asarray(grid.linear_a), initial_linear_a)
    np.testing.assert_allclose(
        np.asarray(grid.native_angular_linear_response),
        initial_response,
    )

    try:
        grid.apply_native_angular_linear_response(mode="unknown")
    except ValueError as exc:
        assert "mode" in str(exc)
    else:
        raise AssertionError("expected coupling-mode validation")
    np.testing.assert_allclose(np.asarray(grid.linear_a), initial_linear_a)
    np.testing.assert_allclose(
        np.asarray(grid.native_angular_linear_response),
        initial_response,
    )


def test_aethergrid_native_angular_linear_feedback_can_drive_linear_velocity_opt_in():
    grid = fdtd.AetherGrid(shape=(3, 3, 3))
    response = np.zeros((3, 3, 3, 3))
    response[1, 1, 1, 2] = 2.0

    grid.apply_native_angular_linear_response(response=response, mode="replace")
    grid.advance_linear_sector()

    np.testing.assert_allclose(
        float(grid.linear_v[1, 1, 1, 2]),
        grid.time_step * 2.0,
    )
    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_A))
    assert not np.any(np.asarray(grid.angular_tau))


def test_aethergrid_step_dynamic_native_torque_feedback_drives_bounded_linear_response():
    def run_dynamic_transport(feedback_mode=None):
        grid = fdtd.AetherGrid(
            shape=(5, 5, 5),
            grid_spacing=8.25e6,
            native_angular_transport=True,
            native_angular_collect_sources=False,
            native_angular_update_clocks=True,
            native_angular_feedback_mode=feedback_mode,
            native_angular_feedback_scale=1.0e-26,
        )
        steps = 32
        grid.ell_t[:, :, :, 0] = 0.1
        grid.ell_p[:, :, :, 0] = 0.1
        grid.angular_inertia_t[:, :, :, 0] = 1.0
        grid.angular_inertia_p[:, :, :, 0] = 1.0
        grid.angular_momentum_t[:, :, :, 0] = 1.0
        grid.angular_momentum_p[:, :, :, 0] = 1.0
        initial_t = np.asarray(grid.angular_momentum_t).copy()
        initial_p = np.asarray(grid.angular_momentum_p).copy()
        initial_rate_t = np.zeros((5, 5, 5, 1))
        initial_rate_p = np.zeros((5, 5, 5, 1))
        initial_rate_t[2, 2, 2, 0] = 0.2
        initial_rate_p[2, 2, 2, 0] = -0.1
        previous_t = initial_t - grid.time_step * initial_rate_t
        previous_p = initial_p - grid.time_step * initial_rate_p
        max_response = 0.0
        max_linear_v = 0.0

        for _ in range(steps):
            current_t = np.asarray(grid.angular_momentum_t).copy()
            current_p = np.asarray(grid.angular_momentum_p).copy()
            grid.update_native_angular_torque(
                previous_momentum_t=previous_t,
                previous_momentum_p=previous_p,
                delta=grid.time_step,
            )
            grid.project_native_angular_torque()

            grid.step()

            response = np.asarray(grid.native_angular_linear_response)
            linear_v = np.asarray(grid.linear_v)
            max_response = max(max_response, float(np.max(np.abs(response))))
            max_linear_v = max(max_linear_v, float(np.max(np.abs(linear_v))))
            assert np.all(np.isfinite(response))
            assert np.all(np.isfinite(linear_v))

            previous_t = current_t
            previous_p = current_p

        return grid, max_response, max_linear_v

    quiet_grid, quiet_response, quiet_linear_v = run_dynamic_transport()
    feedback_grid, feedback_response, feedback_linear_v = run_dynamic_transport(
        feedback_mode="replace"
    )

    np.testing.assert_allclose(
        np.asarray(feedback_grid.angular_momentum_t),
        np.asarray(quiet_grid.angular_momentum_t),
    )
    np.testing.assert_allclose(
        np.asarray(feedback_grid.angular_momentum_p),
        np.asarray(quiet_grid.angular_momentum_p),
    )
    assert quiet_response == 0.0
    assert quiet_linear_v == 0.0
    assert 0.0 < feedback_response < 1.0e-18
    assert 0.0 < feedback_linear_v < 1.0e-5
    assert feedback_grid.time_steps_passed == quiet_grid.time_steps_passed == 32
    assert not np.any(np.asarray(feedback_grid.native_angular_source_t))
    assert not np.any(np.asarray(feedback_grid.native_angular_source_p))


def test_aethergrid_evaluates_linear_charge_flux_candidate():
    grid = fdtd.AetherGrid(shape=(3, 3, 3))
    delta_t = 0.2
    area_measure = np.ones((3, 3, 3, 1)) * 0.25
    grid.linear_v[1, 1, 1] = [2.0, -3.0, 4.0]

    candidate = grid.evaluate_linear_charge_flux_candidate(
        delta_t=delta_t,
        area_measure=area_measure,
        normalization_area=0.5,
        sign=-1.0,
    )

    expected = np.zeros((3, 3, 3, 3))
    expected[1, 1, 1] = (
        -1.0 * delta_t * eta * np.array([2.0, -3.0, 4.0]) * 0.25 / 0.5
    )
    np.testing.assert_allclose(np.asarray(candidate), expected)
    assert candidate is grid.linear_charge_flux_candidate
    assert not np.any(np.asarray(grid.native_angular_charge_flux_t))
    assert not np.any(np.asarray(grid.native_angular_charge_flux_p))
    assert grid.time_steps_passed == 0


def test_aethergrid_linear_charge_flux_requires_reference_area():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))

    try:
        grid.evaluate_linear_charge_flux_candidate(normalization_area=0.0)
    except ValueError as exc:
        assert "normalization_area" in str(exc)
    else:
        raise AssertionError("expected reference-area validation")


def test_aethergrid_evaluates_native_angular_charge_flux_candidates():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))
    delta_t = 0.1
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 5.0
    grid.omega_t[:, :, :, 0] = 3.0
    grid.omega_p[:, :, :, 0] = 7.0
    loop_t = np.ones((2, 2, 2, 1)) * 0.25
    loop_p = np.ones((2, 2, 2, 1)) * 0.5

    q_t, q_p = grid.evaluate_native_angular_charge_flux_candidates(
        delta_t=delta_t,
        loop_measure_t=loop_t,
        loop_measure_p=loop_p,
        angular_normalization=2.0,
        sign_t=-1.0,
        sign_p=1.0,
    )

    np.testing.assert_allclose(
        np.asarray(q_t),
        -1.0 * delta_t * eta * 2.0 * 3.0 * 0.25 / 2.0,
    )
    np.testing.assert_allclose(
        np.asarray(q_p),
        delta_t * eta * 5.0 * 7.0 * 0.5 / 2.0,
    )
    assert q_t is grid.native_angular_charge_flux_t
    assert q_p is grid.native_angular_charge_flux_p
    assert not np.any(np.asarray(grid.linear_charge_flux_candidate))
    assert grid.time_steps_passed == 0


def test_aethergrid_native_angular_charge_flux_requires_normalization():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))

    try:
        grid.evaluate_native_angular_charge_flux_candidates(
            angular_normalization=0.0,
        )
    except ValueError as exc:
        assert "angular_normalization" in str(exc)
    else:
        raise AssertionError("expected angular-normalization validation")


def test_aethergrid_reduces_native_angular_charge_flux_candidates():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))
    grid.native_angular_charge_flux_t[:, :, :, 0] = 4.0
    grid.native_angular_charge_flux_p[:, :, :, 0] = 9.0

    additive = grid.reduce_native_angular_charge_flux(mode="additive")
    geometric = grid.reduce_native_angular_charge_flux(mode="geometric_mean")

    np.testing.assert_allclose(np.asarray(additive), 13.0)
    np.testing.assert_allclose(np.asarray(geometric), 6.0)
    assert geometric is grid.native_angular_charge_reduction
    assert not np.any(np.asarray(grid.angular_A))
    assert not np.any(np.asarray(grid.linear_a))
    assert grid.time_steps_passed == 0


def test_aethergrid_rejects_invalid_native_angular_charge_reduction():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))
    grid.native_angular_charge_flux_t[:, :, :, 0] = -4.0
    grid.native_angular_charge_flux_p[:, :, :, 0] = 9.0

    try:
        grid.reduce_native_angular_charge_flux(mode="geometric_mean")
    except ValueError as exc:
        assert "non-negative" in str(exc)
    else:
        raise AssertionError("expected geometric-mean sign validation")

    try:
        grid.reduce_native_angular_charge_flux(mode="unknown")
    except ValueError as exc:
        assert "mode" in str(exc)
    else:
        raise AssertionError("expected reduction-mode validation")


def test_aethergrid_evaluates_native_angular_transport_residual():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_torque_t[:, :, :, 0] = 0.25
    grid.angular_torque_p[:, :, :, 0] = 0.5
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    source_t = np.ones((4, 4, 4, 1)) * 0.1
    source_p = np.ones((4, 4, 4, 1)) * 0.2

    residual_t, residual_p = grid.evaluate_native_angular_transport_residual(
        source_t=source_t,
        source_p=source_p,
    )
    expected_divergence = np.asarray(div(grid.native_angular_tau))

    np.testing.assert_allclose(
        np.asarray(residual_t),
        0.25 + 2.0 * expected_divergence - source_t,
    )
    np.testing.assert_allclose(
        np.asarray(residual_p),
        0.5 + 3.0 * expected_divergence - source_p,
    )
    assert residual_t is grid.native_angular_transport_residual_t
    assert residual_p is grid.native_angular_transport_residual_p
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_transport_residual_zeroes_when_source_matches_balance():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_torque_t[:, :, :, 0] = 0.25
    grid.angular_torque_p[:, :, :, 0] = 0.5
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0

    metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
    source_t = np.asarray(grid.angular_torque_t) + np.asarray(metric_t)
    source_p = np.asarray(grid.angular_torque_p) + np.asarray(metric_p)

    residual_t, residual_p = grid.evaluate_native_angular_transport_residual(
        source_t=source_t,
        source_p=source_p,
    )

    np.testing.assert_allclose(np.asarray(residual_t), 0.0)
    np.testing.assert_allclose(np.asarray(residual_p), 0.0)


def test_aethergrid_evaluates_native_angular_momentum_transport_rhs():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    source_t = np.ones((4, 4, 4, 1)) * 0.1
    source_p = np.ones((4, 4, 4, 1)) * 0.2

    rhs_t, rhs_p = grid.evaluate_native_angular_momentum_rhs(
        source_t=source_t,
        source_p=source_p,
    )
    expected_divergence = np.asarray(div(grid.native_angular_tau))

    np.testing.assert_allclose(
        np.asarray(rhs_t),
        source_t - 2.0 * expected_divergence,
    )
    np.testing.assert_allclose(
        np.asarray(rhs_p),
        source_p - 3.0 * expected_divergence,
    )
    assert rhs_t is grid.native_angular_momentum_rhs_t
    assert rhs_p is grid.native_angular_momentum_rhs_p
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_predicts_passive_native_angular_momentum_step():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.125
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    source_t = np.ones((4, 4, 4, 1)) * 0.1
    source_p = np.ones((4, 4, 4, 1)) * 0.2
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()

    candidate_t, candidate_p = grid.predict_native_angular_momentum_step(
        delta=delta,
        source_t=source_t,
        source_p=source_p,
    )
    expected_divergence = np.asarray(div(grid.native_angular_tau))
    expected_t = initial_momentum_t + delta * (
        source_t - 2.0 * expected_divergence
    )
    expected_p = initial_momentum_p + delta * (
        source_p - 3.0 * expected_divergence
    )

    np.testing.assert_allclose(np.asarray(candidate_t), expected_t)
    np.testing.assert_allclose(np.asarray(candidate_p), expected_p)
    assert candidate_t is grid.native_angular_momentum_candidate_t
    assert candidate_p is grid.native_angular_momentum_candidate_p
    np.testing.assert_allclose(np.asarray(grid.angular_momentum_t), initial_momentum_t)
    np.testing.assert_allclose(np.asarray(grid.angular_momentum_p), initial_momentum_p)
    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_advances_native_angular_momentum_transport_opt_in():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.125
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    source_t = np.ones((4, 4, 4, 1)) * 0.1
    source_p = np.ones((4, 4, 4, 1)) * 0.2
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()

    momentum_t, momentum_p = grid.advance_native_angular_momentum_transport(
        delta=delta,
        source_t=source_t,
        source_p=source_p,
    )
    expected_divergence = np.asarray(div(grid.native_angular_tau))
    expected_t = initial_momentum_t + delta * (
        source_t - 2.0 * expected_divergence
    )
    expected_p = initial_momentum_p + delta * (
        source_p - 3.0 * expected_divergence
    )

    np.testing.assert_allclose(np.asarray(momentum_t), expected_t)
    np.testing.assert_allclose(np.asarray(momentum_p), expected_p)
    assert momentum_t is grid.angular_momentum_t
    assert momentum_p is grid.angular_momentum_p
    assert momentum_t is grid.native_angular_momentum_candidate_t
    assert momentum_p is grid.native_angular_momentum_candidate_p
    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.angular_A))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_transport_advance_can_collect_source_terms():
    grid = fdtd.AetherGrid(shape=(3, 3, 3))
    delta = 0.2
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    explicit_t = np.ones((3, 3, 3, 1)) * 0.1
    explicit_p = np.ones((3, 3, 3, 1)) * 0.2
    source_t = np.ones((3, 3, 3, 1)) * 0.3
    source_p = np.ones((3, 3, 3, 1)) * 0.4
    boundary_t = np.ones((3, 3, 3, 1)) * 0.5
    boundary_p = np.ones((3, 3, 3, 1)) * 0.6
    grid[1, 1, 1] = NativeAngularExchangeProbe(
        source_t,
        source_p,
        "sources",
    )
    grid[0, :, :] = NativeAngularExchangeProbe(
        boundary_t,
        boundary_p,
        "boundaries",
    )
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()

    momentum_t, momentum_p = grid.advance_native_angular_momentum_transport(
        delta=delta,
        source_t=explicit_t,
        source_p=explicit_p,
        collect_sources=True,
    )

    collected_t = source_t + boundary_t
    collected_p = source_p + boundary_p
    np.testing.assert_allclose(np.asarray(grid.native_angular_source_t), collected_t)
    np.testing.assert_allclose(np.asarray(grid.native_angular_source_p), collected_p)
    np.testing.assert_allclose(
        np.asarray(momentum_t),
        initial_momentum_t + delta * (explicit_t + collected_t),
    )
    np.testing.assert_allclose(
        np.asarray(momentum_p),
        initial_momentum_p + delta * (explicit_p + collected_p),
    )
    assert not np.any(np.asarray(grid.native_angular_tau))
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_transport_advance_can_update_native_clocks():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))
    delta = 0.25
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 1.0
    grid.angular_momentum_p[:, :, :, 0] = 2.0
    source_t = np.ones((2, 2, 2, 1)) * 0.5
    source_p = np.ones((2, 2, 2, 1)) * 1.0

    momentum_t, momentum_p = grid.advance_native_angular_momentum_transport(
        delta=delta,
        source_t=source_t,
        source_p=source_p,
        update_clocks=True,
    )

    np.testing.assert_allclose(np.asarray(momentum_t), 1.0 + delta * 0.5)
    np.testing.assert_allclose(np.asarray(momentum_p), 2.0 + delta * 1.0)
    np.testing.assert_allclose(np.asarray(grid.omega_t), np.asarray(momentum_t) / 2.0)
    np.testing.assert_allclose(np.asarray(grid.omega_p), np.asarray(momentum_p) / 4.0)


def test_aethergrid_transport_advance_requires_inertia_for_clock_update():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))
    grid.angular_inertia_t[:, :, :, 0] = 1.0
    grid.angular_inertia_p[:, :, :, 0] = 0.0
    grid.angular_momentum_t[:, :, :, 0] = 1.0
    grid.angular_momentum_p[:, :, :, 0] = 2.0
    initial_t = np.asarray(grid.angular_momentum_t).copy()
    initial_p = np.asarray(grid.angular_momentum_p).copy()

    try:
        grid.advance_native_angular_momentum_transport(update_clocks=True)
    except ValueError as exc:
        assert "inertia" in str(exc)
    else:
        raise AssertionError("expected positive-inertia validation")
    np.testing.assert_allclose(np.asarray(grid.angular_momentum_t), initial_t)
    np.testing.assert_allclose(np.asarray(grid.angular_momentum_p), initial_p)


def test_aethergrid_residual_and_predictor_ignore_collected_sources_by_default():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.125
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.angular_torque_t[:, :, :, 0] = 0.25
    grid.angular_torque_p[:, :, :, 0] = 0.5
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    boundary_t = np.ones((4, 4, 4, 1)) * 0.75
    boundary_p = np.ones((4, 4, 4, 1)) * 1.25
    grid[0, :, :] = NativeAngularExchangeProbe(
        boundary_t,
        boundary_p,
        "boundaries",
    )
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()

    collected_t, collected_p = grid.collect_native_angular_source_terms()
    collected_t = np.asarray(collected_t).copy()
    collected_p = np.asarray(collected_p).copy()
    expected_divergence = np.asarray(div(grid.native_angular_tau))
    metric_t = 2.0 * expected_divergence
    metric_p = 3.0 * expected_divergence

    residual_t, residual_p = grid.evaluate_native_angular_transport_residual()
    np.testing.assert_allclose(np.asarray(residual_t), 0.25 + metric_t)
    np.testing.assert_allclose(np.asarray(residual_p), 0.5 + metric_p)
    np.testing.assert_allclose(np.asarray(grid.native_angular_source_t), collected_t)
    np.testing.assert_allclose(np.asarray(grid.native_angular_source_p), collected_p)

    candidate_t, candidate_p = grid.predict_native_angular_momentum_step(delta=delta)
    np.testing.assert_allclose(
        np.asarray(candidate_t),
        initial_momentum_t - delta * metric_t,
    )
    np.testing.assert_allclose(
        np.asarray(candidate_p),
        initial_momentum_p - delta * metric_p,
    )
    np.testing.assert_allclose(np.asarray(grid.native_angular_source_t), collected_t)
    np.testing.assert_allclose(np.asarray(grid.native_angular_source_p), collected_p)

    residual_t, residual_p = grid.evaluate_native_angular_transport_residual(
        source_t=collected_t,
        source_p=collected_p,
    )
    np.testing.assert_allclose(np.asarray(residual_t), 0.25 + metric_t - collected_t)
    np.testing.assert_allclose(np.asarray(residual_p), 0.5 + metric_p - collected_p)


def test_aethergrid_passive_momentum_predictor_balanced_64_step_benchmark():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.05
    steps = 64
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
    source_t = np.asarray(metric_t).copy()
    source_p = np.asarray(metric_p).copy()
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()

    for _ in range(steps):
        candidate_t, candidate_p = grid.predict_native_angular_momentum_step(
            delta=delta,
            source_t=source_t,
            source_p=source_p,
        )
        np.testing.assert_allclose(np.asarray(candidate_t), initial_momentum_t)
        np.testing.assert_allclose(np.asarray(candidate_p), initial_momentum_p)
        assert np.all(np.isfinite(np.asarray(candidate_t)))
        assert np.all(np.isfinite(np.asarray(candidate_p)))

        # This assignment simulates a candidate-only update loop inside the
        # benchmark. The production stepper still does not call this predictor.
        grid.angular_momentum_t = candidate_t
        grid.angular_momentum_p = candidate_p

    assert grid.time_steps_passed == 0
    np.testing.assert_allclose(np.asarray(grid.angular_momentum_t), initial_momentum_t)
    np.testing.assert_allclose(np.asarray(grid.angular_momentum_p), initial_momentum_p)
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_native_angular_transport_advance_balanced_64_step_benchmark():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.05
    steps = 64
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
    source_t = np.asarray(metric_t).copy()
    source_p = np.asarray(metric_p).copy()
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()

    for _ in range(steps):
        momentum_t, momentum_p = grid.advance_native_angular_momentum_transport(
            delta=delta,
            source_t=source_t,
            source_p=source_p,
        )
        np.testing.assert_allclose(np.asarray(momentum_t), initial_momentum_t)
        np.testing.assert_allclose(np.asarray(momentum_p), initial_momentum_p)
        assert np.all(np.isfinite(np.asarray(momentum_t)))
        assert np.all(np.isfinite(np.asarray(momentum_p)))

    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.angular_A))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_passive_momentum_predictor_damping_64_step_benchmark():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.05
    steps = 64
    damping_t = 0.2
    damping_p = 0.125
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
    metric_t = np.asarray(metric_t).copy()
    metric_p = np.asarray(metric_p).copy()
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()
    previous_t = initial_momentum_t
    previous_p = initial_momentum_p

    for step in range(steps):
        source_t = metric_t - damping_t * np.asarray(grid.angular_momentum_t)
        source_p = metric_p - damping_p * np.asarray(grid.angular_momentum_p)
        candidate_t, candidate_p = grid.predict_native_angular_momentum_step(
            delta=delta,
            source_t=source_t,
            source_p=source_p,
        )
        candidate_t = np.asarray(candidate_t)
        candidate_p = np.asarray(candidate_p)
        expected_t = initial_momentum_t * (1.0 - damping_t * delta) ** (step + 1)
        expected_p = initial_momentum_p * (1.0 - damping_p * delta) ** (step + 1)

        np.testing.assert_allclose(candidate_t, expected_t)
        np.testing.assert_allclose(candidate_p, expected_p)
        assert np.all(candidate_t <= previous_t)
        assert np.all(candidate_p <= previous_p)
        assert np.all(candidate_t >= 0.0)
        assert np.all(candidate_p >= 0.0)
        assert np.all(np.isfinite(candidate_t))
        assert np.all(np.isfinite(candidate_p))

        # Candidate-only benchmark loop; production dynamics still do not call
        # or apply the predictor.
        grid.angular_momentum_t = grid.native_angular_momentum_candidate_t
        grid.angular_momentum_p = grid.native_angular_momentum_candidate_p
        previous_t = candidate_t
        previous_p = candidate_p

    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_passive_momentum_predictor_sponge_64_step_benchmark():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.05
    steps = 64
    damping_t = np.zeros((4, 4, 4, 1))
    damping_p = np.zeros((4, 4, 4, 1))
    damping_t[0, :, :, 0] = 0.25
    damping_t[-1, :, :, 0] = 0.25
    damping_p[0, :, :, 0] = 0.15
    damping_p[-1, :, :, 0] = 0.15
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
    metric_t = np.asarray(metric_t).copy()
    metric_p = np.asarray(metric_p).copy()
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()

    for step in range(steps):
        source_t = metric_t - damping_t * np.asarray(grid.angular_momentum_t)
        source_p = metric_p - damping_p * np.asarray(grid.angular_momentum_p)
        candidate_t, candidate_p = grid.predict_native_angular_momentum_step(
            delta=delta,
            source_t=source_t,
            source_p=source_p,
        )
        candidate_t = np.asarray(candidate_t)
        candidate_p = np.asarray(candidate_p)
        expected_t = initial_momentum_t * (1.0 - damping_t * delta) ** (step + 1)
        expected_p = initial_momentum_p * (1.0 - damping_p * delta) ** (step + 1)

        np.testing.assert_allclose(candidate_t, expected_t)
        np.testing.assert_allclose(candidate_p, expected_p)
        np.testing.assert_allclose(candidate_t[1:3], initial_momentum_t[1:3])
        np.testing.assert_allclose(candidate_p[1:3], initial_momentum_p[1:3])
        assert np.all(candidate_t[[0, -1]] <= initial_momentum_t[[0, -1]])
        assert np.all(candidate_p[[0, -1]] <= initial_momentum_p[[0, -1]])
        assert np.all(candidate_t >= 0.0)
        assert np.all(candidate_p >= 0.0)
        assert np.all(np.isfinite(candidate_t))
        assert np.all(np.isfinite(candidate_p))

        # Candidate-only sponge accounting benchmark; no production boundary
        # hook or native angular update is invoked.
        grid.angular_momentum_t = grid.native_angular_momentum_candidate_t
        grid.angular_momentum_p = grid.native_angular_momentum_candidate_p

    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aether_angular_sponge_boundary_exposes_exchange_terms():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )

    source_t, source_p = grid.collect_native_angular_source_terms()

    expected_t = np.zeros((4, 4, 4, 1))
    expected_p = np.zeros((4, 4, 4, 1))
    expected_t[0, :, :, 0] = -0.25 * 0.75
    expected_p[0, :, :, 0] = -0.15 * 1.25
    np.testing.assert_allclose(np.asarray(source_t), expected_t)
    np.testing.assert_allclose(np.asarray(source_p), expected_p)
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aether_angular_sponge_boundary_declares_source_accounting_contract():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )

    contract = grid.boundaries[0].native_angular_boundary_contract()

    assert contract["scope"] == "source_accounting"
    assert contract["exchange"] == "local_damping"
    assert contract["returns"] == (
        "native_angular_source_t",
        "native_angular_source_p",
    )
    assert contract["reads"] == ("angular_momentum_t", "angular_momentum_p")
    assert contract["mutates"] == ()
    assert contract["physical_boundary_law"] is False
    assert contract["requires_outer_face"] is True
    assert contract["boundary_axis"] == "x"
    assert contract["boundary_side"] == "low"


def test_aether_angular_no_exchange_boundary_exposes_zero_exchange_terms():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid[0, :, :] = fdtd.AetherAngularNoExchangeBoundary()

    source_t, source_p = grid.collect_native_angular_source_terms()

    assert not np.any(np.asarray(source_t))
    assert not np.any(np.asarray(source_p))
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aether_angular_no_exchange_boundary_declares_source_accounting_contract():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid[0, :, :] = fdtd.AetherAngularNoExchangeBoundary()

    contract = grid.boundaries[0].native_angular_boundary_contract()

    assert contract["scope"] == "source_accounting"
    assert contract["exchange"] == "none"
    assert contract["returns"] == (
        "native_angular_source_t",
        "native_angular_source_p",
    )
    assert contract["reads"] == ()
    assert contract["mutates"] == ()
    assert contract["physical_boundary_law"] is False
    assert contract["requires_outer_face"] is True
    assert contract["boundary_axis"] == "x"
    assert contract["boundary_side"] == "low"


def test_aether_angular_no_exchange_boundary_has_zero_exchange_power():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid[0, :, :] = fdtd.AetherAngularNoExchangeBoundary()

    source_t, source_p = grid.collect_native_angular_source_terms()
    power_t, power_p, total_power = grid.evaluate_native_angular_exchange_power(
        source_t=source_t,
        source_p=source_p,
    )

    assert not np.any(np.asarray(power_t))
    assert not np.any(np.asarray(power_p))
    assert float(total_power) == 0.0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_boundary_source_collection_mutates_only_source_buffers():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.omega_t[:, :, :, 0] = 0.4
    grid.omega_p[:, :, :, 0] = 0.6
    grid.theta_t[:, :, :, 0] = 0.04
    grid.theta_p[:, :, :, 0] = 0.06
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 3.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.angular_momentum_t[-1, :, :, 0] = 0.5
    grid.angular_momentum_p[-1, :, :, 0] = 0.25
    grid.angular_torque_t[:, :, :, 0] = 0.2
    grid.angular_torque_p[:, :, :, 0] = 0.3
    grid.native_angular_tau[:, :, :, 0] = 0.1
    grid.native_angular_transport_residual_t[:, :, :, 0] = 0.7
    grid.native_angular_transport_residual_p[:, :, :, 0] = 0.9
    grid.native_angular_momentum_candidate_t[:, :, :, 0] = 1.1
    grid.native_angular_momentum_candidate_p[:, :, :, 0] = 1.3
    snapshots = {
        "omega_t": np.asarray(grid.omega_t).copy(),
        "omega_p": np.asarray(grid.omega_p).copy(),
        "theta_t": np.asarray(grid.theta_t).copy(),
        "theta_p": np.asarray(grid.theta_p).copy(),
        "angular_inertia_t": np.asarray(grid.angular_inertia_t).copy(),
        "angular_inertia_p": np.asarray(grid.angular_inertia_p).copy(),
        "angular_momentum_t": np.asarray(grid.angular_momentum_t).copy(),
        "angular_momentum_p": np.asarray(grid.angular_momentum_p).copy(),
        "angular_torque_t": np.asarray(grid.angular_torque_t).copy(),
        "angular_torque_p": np.asarray(grid.angular_torque_p).copy(),
        "native_angular_tau": np.asarray(grid.native_angular_tau).copy(),
        "native_angular_transport_residual_t": np.asarray(
            grid.native_angular_transport_residual_t
        ).copy(),
        "native_angular_transport_residual_p": np.asarray(
            grid.native_angular_transport_residual_p
        ).copy(),
        "native_angular_momentum_candidate_t": np.asarray(
            grid.native_angular_momentum_candidate_t
        ).copy(),
        "native_angular_momentum_candidate_p": np.asarray(
            grid.native_angular_momentum_candidate_p
        ).copy(),
        "linear_a": np.asarray(grid.linear_a).copy(),
        "angular_tau": np.asarray(grid.angular_tau).copy(),
    }
    grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )
    grid[-1, :, :] = fdtd.AetherAngularReflectiveBoundary(
        response_rate_t=2.0,
        response_rate_p=3.0,
        sign_t=-1.0,
        sign_p=1.0,
    )

    source_t, source_p = grid.collect_native_angular_source_terms()

    expected_t = np.zeros((4, 4, 4, 1))
    expected_p = np.zeros((4, 4, 4, 1))
    expected_t[0, :, :, 0] = -0.25 * 0.75
    expected_p[0, :, :, 0] = -0.15 * 1.25
    expected_t[-1, :, :, 0] = 2.0 * (-0.75 - 0.5)
    expected_p[-1, :, :, 0] = 3.0 * (1.25 - 0.25)
    np.testing.assert_allclose(np.asarray(source_t), expected_t)
    np.testing.assert_allclose(np.asarray(source_p), expected_p)

    for name, snapshot in snapshots.items():
        np.testing.assert_allclose(np.asarray(getattr(grid, name)), snapshot)
    assert grid.time_steps_passed == 0


def test_aethergrid_passive_momentum_predictor_no_exchange_boundary_hook_64_step():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.05
    steps = 64
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
    metric_t = np.asarray(metric_t).copy()
    metric_p = np.asarray(metric_p).copy()
    grid[1, 1, 1] = NativeAngularExchangeProbe(metric_t, metric_p, "sources")
    grid[0, :, :] = fdtd.AetherAngularNoExchangeBoundary()
    grid[-1, :, :] = fdtd.AetherAngularNoExchangeBoundary()
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()

    for _ in range(steps):
        source_t, source_p = grid.collect_native_angular_source_terms()
        candidate_t, candidate_p = grid.predict_native_angular_momentum_step(
            delta=delta,
            source_t=source_t,
            source_p=source_p,
        )
        np.testing.assert_allclose(np.asarray(candidate_t), initial_momentum_t)
        np.testing.assert_allclose(np.asarray(candidate_p), initial_momentum_p)
        assert np.all(np.isfinite(np.asarray(candidate_t)))
        assert np.all(np.isfinite(np.asarray(candidate_p)))

        # Candidate-only boundary-hook accounting benchmark; production
        # dynamics still do not call or apply the predictor.
        grid.angular_momentum_t = grid.native_angular_momentum_candidate_t
        grid.angular_momentum_p = grid.native_angular_momentum_candidate_p

    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_transport_advance_no_exchange_boundary_hook_64_step():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.05
    steps = 64
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
    metric_t = np.asarray(metric_t).copy()
    metric_p = np.asarray(metric_p).copy()
    grid[1, 1, 1] = NativeAngularExchangeProbe(metric_t, metric_p, "sources")
    grid[0, :, :] = fdtd.AetherAngularNoExchangeBoundary()
    grid[-1, :, :] = fdtd.AetherAngularNoExchangeBoundary()
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()

    for _ in range(steps):
        momentum_t, momentum_p = grid.advance_native_angular_momentum_transport(
            delta=delta,
            collect_sources=True,
        )
        np.testing.assert_allclose(np.asarray(momentum_t), initial_momentum_t)
        np.testing.assert_allclose(np.asarray(momentum_p), initial_momentum_p)
        assert np.all(np.isfinite(np.asarray(momentum_t)))
        assert np.all(np.isfinite(np.asarray(momentum_p)))

    np.testing.assert_allclose(np.asarray(grid.native_angular_source_t), metric_t)
    np.testing.assert_allclose(np.asarray(grid.native_angular_source_p), metric_p)
    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.angular_A))
    assert not np.any(np.asarray(grid.linear_a))


def test_aether_angular_reflective_boundary_exposes_mirror_exchange_terms():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    response_rate_t = 20.0
    response_rate_p = 10.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.angular_momentum_t[0, :, :, 0] = 0.5
    grid.angular_momentum_p[0, :, :, 0] = 0.25
    grid[0, :, :] = fdtd.AetherAngularReflectiveBoundary(
        response_rate_t=response_rate_t,
        response_rate_p=response_rate_p,
        sign_t=-1.0,
        sign_p=1.0,
    )

    source_t, source_p = grid.collect_native_angular_source_terms()

    expected_t = np.zeros((4, 4, 4, 1))
    expected_p = np.zeros((4, 4, 4, 1))
    expected_t[0, :, :, 0] = response_rate_t * (-0.75 - 0.5)
    expected_p[0, :, :, 0] = response_rate_p * (1.25 - 0.25)
    np.testing.assert_allclose(np.asarray(source_t), expected_t)
    np.testing.assert_allclose(np.asarray(source_p), expected_p)
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aether_angular_reflective_boundary_declares_source_accounting_contract():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid[-1, :, :] = fdtd.AetherAngularReflectiveBoundary(
        response_rate_t=2.0,
        response_rate_p=3.0,
        sign_t=-1.0,
        sign_p=1.0,
    )

    contract = grid.boundaries[0].native_angular_boundary_contract()

    assert contract["scope"] == "source_accounting"
    assert contract["exchange"] == "mirror_relaxation"
    assert contract["returns"] == (
        "native_angular_source_t",
        "native_angular_source_p",
    )
    assert contract["reads"] == ("angular_momentum_t", "angular_momentum_p")
    assert contract["mutates"] == ()
    assert contract["physical_boundary_law"] is False
    assert contract["requires_outer_face"] is True
    assert contract["boundary_axis"] == "x"
    assert contract["boundary_side"] == "high"


def test_aethergrid_collects_native_angular_boundary_contracts():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
        name="sponge_low",
    )
    grid[-1, :, :] = fdtd.AetherAngularReflectiveBoundary(
        response_rate_t=2.0,
        response_rate_p=3.0,
        sign_t=-1.0,
        sign_p=1.0,
        name="mirror_high",
    )

    contracts = grid.collect_native_angular_boundary_contracts()

    assert len(contracts) == 2
    assert contracts[0]["boundary_type"] == "AetherAngularSpongeBoundary"
    assert contracts[0]["boundary_name"] == "sponge_low"
    assert contracts[0]["scope"] == "source_accounting"
    assert contracts[0]["boundary_axis"] == "x"
    assert contracts[0]["boundary_side"] == "low"
    assert contracts[0]["physical_boundary_law"] is False
    assert contracts[1]["boundary_type"] == "AetherAngularReflectiveBoundary"
    assert contracts[1]["boundary_name"] == "mirror_high"
    assert contracts[1]["exchange"] == "mirror_relaxation"
    assert contracts[1]["boundary_axis"] == "x"
    assert contracts[1]["boundary_side"] == "high"
    assert contracts[1]["physical_boundary_law"] is False
    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.native_angular_source_t))
    assert not np.any(np.asarray(grid.native_angular_source_p))


def test_aethergrid_validates_current_native_angular_boundary_contracts():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
        name="sponge_low",
    )
    grid[-1, :, :] = fdtd.AetherAngularReflectiveBoundary(
        response_rate_t=2.0,
        response_rate_p=3.0,
        sign_t=-1.0,
        sign_p=1.0,
        name="mirror_high",
    )

    contracts = grid.validate_native_angular_boundary_contracts()

    assert len(contracts) == 2
    assert all(contract["scope"] == "source_accounting" for contract in contracts)
    assert not any(contract["physical_boundary_law"] for contract in contracts)
    assert not np.any(np.asarray(grid.native_angular_source_t))
    assert not np.any(np.asarray(grid.native_angular_source_p))


def test_aethergrid_rejects_source_accounting_boundaries_as_physical_laws():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )

    try:
        grid.validate_native_angular_boundary_contracts(require_physical_laws=True)
    except ValueError as exc:
        assert "physical boundary law required" in str(exc)
    else:
        raise AssertionError("expected physical boundary law validation")


def test_aethergrid_rejects_incomplete_physical_boundary_contract():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid[0, :, :] = NativeAngularBoundaryContractProbe(
        {
            "scope": "physical_boundary",
            "returns": ("native_angular_source_t", "native_angular_source_p"),
            "reads": ("angular_momentum_t", "angular_momentum_p"),
            "mutates": ("angular_momentum_t", "angular_momentum_p"),
            "physical_boundary_law": True,
        }
    )

    try:
        grid.validate_native_angular_boundary_contracts(require_physical_laws=True)
    except ValueError as exc:
        message = str(exc)
        assert "physical native angular boundary contract" in message
        assert "boundary_slots" in message
        assert "flux_split" in message
        assert "metric_frame" in message
        assert "energy_balance" in message
        assert "acceptance_test" in message
    else:
        raise AssertionError("expected incomplete physical contract validation")


def test_aether_angular_matched_flux_boundary_candidate_declares_physical_contract():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid[0, :, :] = fdtd.AetherAngularMatchedFluxBoundary(
        absorption_t=0.5,
        absorption_p=0.25,
        name="matched_low",
    )

    contracts = grid.validate_native_angular_boundary_contracts(
        require_physical_laws=True,
    )
    contract = contracts[0]

    assert len(contracts) == 1
    assert contract["boundary_type"] == "AetherAngularMatchedFluxBoundary"
    assert contract["boundary_name"] == "matched_low"
    assert contract["scope"] == "physical_boundary_candidate"
    assert contract["exchange"] == "matched_incident_flux_damping"
    assert contract["physical_boundary_law"] is True
    assert contract["requires_outer_face"] is True
    assert contract["boundary_axis"] == "x"
    assert contract["boundary_side"] == "low"
    assert contract["boundary_slots"] == ("outer_face_cells",)
    assert contract["flux_split"] == (
        "projected_normal_incident_outgoing_channel_weighted"
    )
    assert "angular_e_t/p" in contract["metric_frame"]
    assert contract["falsification_comparator"] == "direct_native_channel_flux"
    assert contract["zero_incident_response"] == "zero_source_exchange"
    assert "nonpositive" in contract["energy_balance"]
    assert (
        contract["acceptance_test"]
        == "test_aether_angular_matched_flux_boundary_candidate_is_passive"
    )


def test_aether_angular_matched_flux_boundary_candidate_is_passive():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 3.0
    grid.angular_momentum_p[:, :, :, 0] = 5.0
    grid.native_angular_tau[0, 0, 0, 0] = 2.0
    grid.native_angular_tau[0, 1, 0, 0] = -3.0
    grid.native_angular_tau[-1, 0, 0, 0] = 5.0
    grid.native_angular_tau[-1, 1, 0, 0] = -7.0
    grid[0, :, :] = fdtd.AetherAngularMatchedFluxBoundary(
        absorption_t=0.5,
        absorption_p=0.25,
        name="matched_low",
    )
    grid[-1, :, :] = fdtd.AetherAngularMatchedFluxBoundary(
        absorption_t=0.5,
        absorption_p=0.25,
        name="matched_high",
    )

    source_t, source_p = grid.collect_native_angular_source_terms()
    _, _, exchange_power = grid.evaluate_native_angular_exchange_power(
        source_t=source_t,
        source_p=source_p,
    )

    expected_t = np.zeros((4, 4, 4, 1))
    expected_p = np.zeros((4, 4, 4, 1))
    expected_t[0, 0, 0, 0] = -0.5 * 2.0
    expected_t[-1, 1, 0, 0] = -0.5 * 7.0
    np.testing.assert_allclose(np.asarray(source_t), expected_t)
    np.testing.assert_allclose(np.asarray(source_p), expected_p)
    assert float(exchange_power) <= 0.0
    assert not np.any(np.asarray(grid.linear_a))
    assert not np.any(np.asarray(grid.angular_A))


def test_aether_angular_matched_flux_boundary_candidate_no_incident_no_exchange():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 3.0
    grid.angular_momentum_p[:, :, :, 0] = 5.0
    grid.native_angular_tau[0, 0, 0, 0] = -2.0
    grid.native_angular_tau[0, 1, 0, 0] = -3.0
    grid.native_angular_tau[-1, 0, 0, 0] = 5.0
    grid.native_angular_tau[-1, 1, 0, 0] = 7.0
    grid[0, :, :] = fdtd.AetherAngularMatchedFluxBoundary(
        absorption_t=0.5,
        absorption_p=0.25,
        name="matched_low",
    )
    grid[-1, :, :] = fdtd.AetherAngularMatchedFluxBoundary(
        absorption_t=0.5,
        absorption_p=0.25,
        name="matched_high",
    )

    fluxes = grid.collect_native_angular_boundary_channel_fluxes(
        contracts=grid.validate_native_angular_boundary_contracts(
            require_physical_laws=True,
        ),
    )
    source_t, source_p = grid.collect_native_angular_source_terms()
    _, _, exchange_power = grid.evaluate_native_angular_exchange_power(
        source_t=source_t,
        source_p=source_p,
    )

    assert sum(float(flux["incident_total"]) for flux in fluxes) == 0.0
    assert sum(float(flux["outgoing_total"]) for flux in fluxes) > 0.0
    assert not np.any(np.asarray(source_t))
    assert not np.any(np.asarray(source_p))
    assert float(exchange_power) == 0.0
    assert not np.any(np.asarray(grid.linear_a))
    assert not np.any(np.asarray(grid.angular_A))


def test_aether_angular_direct_matched_flux_boundary_candidate_declares_physical_contract():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid[0, :, :] = fdtd.AetherAngularDirectMatchedFluxBoundary(
        absorption_t=0.5,
        absorption_p=0.25,
        name="direct_low",
    )

    contracts = grid.validate_native_angular_boundary_contracts(
        require_physical_laws=True,
    )
    contract = contracts[0]

    assert len(contracts) == 1
    assert contract["boundary_type"] == "AetherAngularDirectMatchedFluxBoundary"
    assert contract["boundary_name"] == "direct_low"
    assert contract["scope"] == "physical_boundary_candidate"
    assert contract["exchange"] == "direct_channel_matched_flux_damping"
    assert contract["physical_boundary_law"] is True
    assert contract["boundary_axis"] == "x"
    assert contract["boundary_side"] == "low"
    assert contract["boundary_slots"] == ("outer_face_cells",)
    assert contract["flux_split"] == "direct_native_channel_incident_outgoing"
    assert "angular_e_t/p" in contract["metric_frame"]
    assert contract["compares_with"] == "matched_incident_flux_damping"
    assert contract["zero_incident_response"] == "zero_source_exchange"
    assert "nonpositive" in contract["energy_balance"]
    assert "direct_matched_flux" in contract["acceptance_test"]


def test_aether_angular_direct_matched_flux_boundary_candidate_is_passive():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 3.0
    grid.angular_momentum_p[:, :, :, 0] = 5.0
    grid.angular_torque_t[0, 0, 0, 0] = 2.0
    grid.angular_torque_t[0, 1, 0, 0] = -3.0
    grid.angular_torque_t[-1, 0, 0, 0] = 5.0
    grid.angular_torque_t[-1, 1, 0, 0] = -7.0
    grid[0, :, :] = fdtd.AetherAngularDirectMatchedFluxBoundary(
        absorption_t=0.5,
        absorption_p=0.25,
        name="direct_low",
    )
    grid[-1, :, :] = fdtd.AetherAngularDirectMatchedFluxBoundary(
        absorption_t=0.5,
        absorption_p=0.25,
        name="direct_high",
    )

    source_t, source_p = grid.collect_native_angular_source_terms()
    _, _, exchange_power = grid.evaluate_native_angular_exchange_power(
        source_t=source_t,
        source_p=source_p,
    )

    expected_t = np.zeros((4, 4, 4, 1))
    expected_p = np.zeros((4, 4, 4, 1))
    expected_t[0, 0, 0, 0] = -0.5 * 2.0
    expected_t[-1, 1, 0, 0] = -0.5 * 7.0
    np.testing.assert_allclose(np.asarray(source_t), expected_t)
    np.testing.assert_allclose(np.asarray(source_p), expected_p)
    assert float(exchange_power) <= 0.0
    assert not np.any(np.asarray(grid.linear_a))
    assert not np.any(np.asarray(grid.angular_A))


def test_aether_angular_direct_matched_flux_boundary_candidate_no_incident_no_exchange():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 3.0
    grid.angular_momentum_p[:, :, :, 0] = 5.0
    grid.angular_torque_t[0, 0, 0, 0] = -2.0
    grid.angular_torque_t[0, 1, 0, 0] = -3.0
    grid.angular_torque_t[-1, 0, 0, 0] = 5.0
    grid.angular_torque_t[-1, 1, 0, 0] = 7.0
    grid[0, :, :] = fdtd.AetherAngularDirectMatchedFluxBoundary(
        absorption_t=0.5,
        absorption_p=0.25,
        name="direct_low",
    )
    grid[-1, :, :] = fdtd.AetherAngularDirectMatchedFluxBoundary(
        absorption_t=0.5,
        absorption_p=0.25,
        name="direct_high",
    )

    source_t, source_p = grid.collect_native_angular_source_terms()
    _, _, exchange_power = grid.evaluate_native_angular_exchange_power(
        source_t=source_t,
        source_p=source_p,
    )

    assert not np.any(np.asarray(source_t))
    assert not np.any(np.asarray(source_p))
    assert float(exchange_power) == 0.0


def test_aether_angular_direct_matched_flux_boundary_candidate_absorbs_counterpropagating_channel_flux():
    matched_grid = fdtd.AetherGrid(shape=(4, 4, 4))
    direct_grid = fdtd.AetherGrid(shape=(4, 4, 4))
    inv_sqrt_2 = 1.0 / np.sqrt(2.0)
    for grid in (matched_grid, direct_grid):
        grid.angular_inertia_t[:, :, :, 0] = 1.0
        grid.angular_inertia_p[:, :, :, 0] = 1.0
        grid.angular_momentum_t[:, :, :, 0] = 1.0
        grid.angular_momentum_p[:, :, :, 0] = 1.0
        grid.configure_native_angular_frame(
            e_t=[inv_sqrt_2, inv_sqrt_2, 0.0],
            e_p=[inv_sqrt_2, -inv_sqrt_2, 0.0],
        )
        grid.angular_torque_t[0, 0, 0, 0] = -2.0 * np.sqrt(2.0)
        grid.angular_torque_p[0, 0, 0, 0] = 2.0 * np.sqrt(2.0)
        grid.project_native_angular_torque()

    matched_grid[0, :, :] = fdtd.AetherAngularMatchedFluxBoundary(
        absorption_t=0.5,
        absorption_p=0.5,
    )
    direct_grid[0, :, :] = fdtd.AetherAngularDirectMatchedFluxBoundary(
        absorption_t=0.5,
        absorption_p=0.5,
    )

    matched_t, matched_p = matched_grid.collect_native_angular_source_terms()
    direct_t, direct_p = direct_grid.collect_native_angular_source_terms()
    _, _, matched_power = matched_grid.evaluate_native_angular_exchange_power(
        source_t=matched_t,
        source_p=matched_p,
    )
    _, _, direct_power = direct_grid.evaluate_native_angular_exchange_power(
        source_t=direct_t,
        source_p=direct_p,
    )

    assert not np.any(np.asarray(matched_t))
    assert not np.any(np.asarray(matched_p))
    assert float(matched_power) == 0.0
    np.testing.assert_allclose(float(direct_t[0, 0, 0, 0]), 0.0)
    np.testing.assert_allclose(float(direct_p[0, 0, 0, 0]), -1.0)
    assert float(direct_power) < 0.0
    assert not np.any(np.asarray(direct_grid.linear_a))
    assert not np.any(np.asarray(direct_grid.angular_A))


def test_aether_angular_matched_flux_boundary_candidate_uses_channel_frame():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 3.0
    grid.angular_momentum_p[:, :, :, 0] = 5.0
    grid.native_angular_tau[0, -1, 0, 1] = 4.0
    grid.native_angular_tau[1, -1, 0, 1] = -6.0
    grid[:, -1, :] = fdtd.AetherAngularMatchedFluxBoundary(
        absorption_t=0.5,
        absorption_p=0.25,
        name="matched_y_high",
    )

    source_t, source_p = grid.collect_native_angular_source_terms()
    _, _, exchange_power = grid.evaluate_native_angular_exchange_power(
        source_t=source_t,
        source_p=source_p,
    )

    expected_t = np.zeros((4, 4, 4, 1))
    expected_p = np.zeros((4, 4, 4, 1))
    expected_p[1, -1, 0, 0] = -0.25 * 6.0
    np.testing.assert_allclose(np.asarray(source_t), expected_t)
    np.testing.assert_allclose(np.asarray(source_p), expected_p)
    assert float(exchange_power) <= 0.0


def test_aethergrid_evaluates_native_angular_boundary_flux_split():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.native_angular_tau[0, 0, 0, 0] = 2.0
    grid.native_angular_tau[0, 1, 0, 0] = -3.0
    grid.native_angular_tau[-1, 0, 0, 0] = 5.0
    grid.native_angular_tau[-1, 1, 0, 0] = -7.0

    low = grid.evaluate_native_angular_boundary_flux(axis="x", side="low")

    assert low["axis"] == "x"
    assert low["side"] == "low"
    assert low["normal_sign"] == -1.0
    np.testing.assert_allclose(float(low["incident_total"]), 2.0)
    np.testing.assert_allclose(float(low["outgoing_total"]), 3.0)
    np.testing.assert_allclose(float(low["net_flux"]), 1.0)
    np.testing.assert_allclose(
        float(grid.native_angular_boundary_incident_flux[0, 0, 0, 0]),
        2.0,
    )
    np.testing.assert_allclose(
        float(grid.native_angular_boundary_outgoing_flux[0, 1, 0, 0]),
        3.0,
    )

    high = grid.evaluate_native_angular_boundary_flux(axis=0, side="high")

    assert high["axis"] == "x"
    assert high["side"] == "high"
    assert high["normal_sign"] == 1.0
    np.testing.assert_allclose(float(high["incident_total"]), 7.0)
    np.testing.assert_allclose(float(high["outgoing_total"]), 5.0)
    np.testing.assert_allclose(float(high["net_flux"]), -2.0)
    np.testing.assert_allclose(
        float(grid.native_angular_boundary_incident_flux[-1, 1, 0, 0]),
        7.0,
    )
    np.testing.assert_allclose(
        float(grid.native_angular_boundary_outgoing_flux[-1, 0, 0, 0]),
        5.0,
    )


def test_aethergrid_collects_native_angular_boundary_fluxes_from_contracts():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.native_angular_tau[0, 0, 0, 0] = 2.0
    grid.native_angular_tau[0, 1, 0, 0] = -3.0
    grid.native_angular_tau[-1, 0, 0, 0] = 5.0
    grid.native_angular_tau[-1, 1, 0, 0] = -7.0
    grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
        name="sponge_low",
    )
    grid[-1, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
        name="sponge_high",
    )

    fluxes = grid.collect_native_angular_boundary_fluxes()

    assert len(fluxes) == 2
    assert fluxes[0]["boundary_name"] == "sponge_low"
    assert fluxes[0]["axis"] == "x"
    assert fluxes[0]["side"] == "low"
    np.testing.assert_allclose(float(fluxes[0]["incident_total"]), 2.0)
    np.testing.assert_allclose(float(fluxes[0]["outgoing_total"]), 3.0)
    np.testing.assert_allclose(float(fluxes[0]["net_flux"]), 1.0)
    assert fluxes[1]["boundary_name"] == "sponge_high"
    assert fluxes[1]["axis"] == "x"
    assert fluxes[1]["side"] == "high"
    np.testing.assert_allclose(float(fluxes[1]["incident_total"]), 7.0)
    np.testing.assert_allclose(float(fluxes[1]["outgoing_total"]), 5.0)
    np.testing.assert_allclose(float(fluxes[1]["net_flux"]), -2.0)
    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.native_angular_source_t))
    assert not np.any(np.asarray(grid.native_angular_source_p))


def test_aethergrid_evaluates_native_angular_boundary_channel_flux_split():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.native_angular_tau[0, 0, 0, 0] = 2.0
    grid.native_angular_tau[0, 1, 0, 0] = -3.0
    grid.native_angular_tau[0, -1, 0, 1] = 4.0
    grid.native_angular_tau[1, -1, 0, 1] = -6.0

    x_low = grid.evaluate_native_angular_boundary_channel_flux(
        axis="x",
        side="low",
    )

    np.testing.assert_allclose(float(x_low["incident_total_t"]), 2.0)
    np.testing.assert_allclose(float(x_low["incident_total_p"]), 0.0)
    np.testing.assert_allclose(float(x_low["outgoing_total_t"]), 3.0)
    np.testing.assert_allclose(float(x_low["outgoing_total_p"]), 0.0)
    np.testing.assert_allclose(
        float(grid.native_angular_boundary_incident_flux_t[0, 0, 0, 0]),
        2.0,
    )
    np.testing.assert_allclose(
        float(grid.native_angular_boundary_outgoing_flux_p[0, 1, 0, 0]),
        0.0,
    )

    y_high = grid.evaluate_native_angular_boundary_channel_flux(
        axis="y",
        side="high",
    )

    np.testing.assert_allclose(float(y_high["incident_total_t"]), 0.0)
    np.testing.assert_allclose(float(y_high["incident_total_p"]), 6.0)
    np.testing.assert_allclose(float(y_high["outgoing_total_t"]), 0.0)
    np.testing.assert_allclose(float(y_high["outgoing_total_p"]), 4.0)
    np.testing.assert_allclose(
        float(grid.native_angular_boundary_incident_flux_p[1, -1, 0, 0]),
        6.0,
    )
    np.testing.assert_allclose(
        float(grid.native_angular_boundary_outgoing_flux_t[0, -1, 0, 0]),
        0.0,
    )


def test_aethergrid_collects_native_angular_boundary_channel_fluxes_from_contracts():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.native_angular_tau[0, 0, 0, 0] = 2.0
    grid.native_angular_tau[0, 1, 0, 0] = -3.0
    grid.native_angular_tau[0, -1, 0, 1] = 4.0
    grid.native_angular_tau[1, -1, 0, 1] = -6.0
    grid[0, :, :] = fdtd.AetherAngularNoExchangeBoundary(name="x_low")
    grid[:, -1, :] = fdtd.AetherAngularNoExchangeBoundary(name="y_high")

    fluxes = grid.collect_native_angular_boundary_channel_fluxes()

    assert len(fluxes) == 2
    assert fluxes[0]["boundary_name"] == "x_low"
    np.testing.assert_allclose(float(fluxes[0]["incident_total_t"]), 2.0)
    np.testing.assert_allclose(float(fluxes[0]["incident_total_p"]), 0.0)
    np.testing.assert_allclose(float(fluxes[0]["outgoing_total_t"]), 3.0)
    np.testing.assert_allclose(float(fluxes[0]["outgoing_total_p"]), 0.0)
    assert fluxes[1]["boundary_name"] == "y_high"
    np.testing.assert_allclose(float(fluxes[1]["incident_total_t"]), 0.0)
    np.testing.assert_allclose(float(fluxes[1]["incident_total_p"]), 6.0)
    np.testing.assert_allclose(float(fluxes[1]["outgoing_total_t"]), 0.0)
    np.testing.assert_allclose(float(fluxes[1]["outgoing_total_p"]), 4.0)
    assert not np.any(np.asarray(grid.native_angular_source_t))
    assert not np.any(np.asarray(grid.native_angular_source_p))


def test_aethergrid_direct_native_channel_flux_matches_axis_aligned_projection():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.angular_torque_t[0, 0, 0, 0] = -2.0
    grid.angular_torque_t[0, 1, 0, 0] = 3.0
    grid.project_native_angular_torque()

    weighted = grid.evaluate_native_angular_boundary_channel_flux(
        axis="x",
        side="low",
    )
    direct = grid.evaluate_native_angular_boundary_direct_channel_flux(
        axis="x",
        side="low",
    )

    np.testing.assert_allclose(
        float(direct["incident_total_t"]),
        float(weighted["incident_total_t"]),
    )
    np.testing.assert_allclose(
        float(direct["outgoing_total_t"]),
        float(weighted["outgoing_total_t"]),
    )
    np.testing.assert_allclose(float(direct["incident_total_p"]), 0.0)
    np.testing.assert_allclose(float(direct["outgoing_total_p"]), 0.0)
    np.testing.assert_allclose(
        float(grid.native_angular_boundary_direct_incident_flux_t[0, 1, 0, 0]),
        3.0,
    )
    np.testing.assert_allclose(
        float(grid.native_angular_boundary_direct_outgoing_flux_t[0, 0, 0, 0]),
        2.0,
    )


def test_aethergrid_direct_native_channel_flux_exposes_projected_cancellation():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    inv_sqrt_2 = 1.0 / np.sqrt(2.0)
    grid.configure_native_angular_frame(
        e_t=[inv_sqrt_2, inv_sqrt_2, 0.0],
        e_p=[inv_sqrt_2, -inv_sqrt_2, 0.0],
    )
    grid.angular_torque_t[0, 0, 0, 0] = -2.0 * np.sqrt(2.0)
    grid.angular_torque_p[0, 0, 0, 0] = 2.0 * np.sqrt(2.0)
    grid.project_native_angular_torque()

    weighted = grid.evaluate_native_angular_boundary_channel_flux(
        axis="x",
        side="low",
    )
    direct = grid.evaluate_native_angular_boundary_direct_channel_flux(
        axis="x",
        side="low",
    )

    np.testing.assert_allclose(float(weighted["incident_total"]), 0.0, atol=1e-12)
    np.testing.assert_allclose(float(weighted["outgoing_total"]), 0.0, atol=1e-12)
    np.testing.assert_allclose(float(direct["net_flux"]), 0.0, atol=1e-12)
    np.testing.assert_allclose(float(direct["outgoing_total_t"]), 2.0)
    np.testing.assert_allclose(float(direct["incident_total_p"]), 2.0)
    np.testing.assert_allclose(float(direct["incident_total"]), 2.0)
    np.testing.assert_allclose(float(direct["outgoing_total"]), 2.0)
    np.testing.assert_allclose(
        float(grid.native_angular_boundary_direct_outgoing_flux_t[0, 0, 0, 0]),
        2.0,
    )
    np.testing.assert_allclose(
        float(grid.native_angular_boundary_direct_incident_flux_p[0, 0, 0, 0]),
        2.0,
    )


def test_aethergrid_collects_direct_native_channel_fluxes_from_contracts():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.angular_torque_t[0, 0, 0, 0] = -2.0
    grid.angular_torque_p[1, -1, 0, 0] = -6.0
    grid[0, :, :] = fdtd.AetherAngularNoExchangeBoundary(name="x_low")
    grid[:, -1, :] = fdtd.AetherAngularNoExchangeBoundary(name="y_high")

    fluxes = grid.collect_native_angular_boundary_direct_channel_fluxes()

    assert len(fluxes) == 2
    assert fluxes[0]["boundary_name"] == "x_low"
    np.testing.assert_allclose(float(fluxes[0]["outgoing_total_t"]), 2.0)
    np.testing.assert_allclose(float(fluxes[0]["incident_total_p"]), 0.0)
    assert fluxes[1]["boundary_name"] == "y_high"
    np.testing.assert_allclose(float(fluxes[1]["incident_total_p"]), 6.0)
    np.testing.assert_allclose(float(fluxes[1]["outgoing_total_t"]), 0.0)
    assert not np.any(np.asarray(grid.native_angular_source_t))
    assert not np.any(np.asarray(grid.native_angular_source_p))


def test_aethergrid_native_angular_boundary_flux_validates_inputs():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))

    for axis, side, message in (
        ("q", "low", "axis must be"),
        ("x", "middle", "side must be"),
    ):
        try:
            grid.evaluate_native_angular_boundary_flux(axis=axis, side=side)
        except ValueError as exc:
            assert message in str(exc)
        else:
            raise AssertionError("expected boundary flux validation")

    try:
        grid.evaluate_native_angular_boundary_flux(
            axis="x",
            side="low",
            field=np.zeros((4, 4, 4, 1)),
        )
    except ValueError as exc:
        assert "boundary flux field must have shape" in str(exc)
    else:
        raise AssertionError("expected boundary flux field validation")

    try:
        grid.evaluate_native_angular_boundary_direct_channel_flux(
            axis="x",
            side="low",
            torque_t=np.zeros((4, 4, 4, 3)),
        )
    except ValueError as exc:
        assert "torque_t" in str(exc)
    else:
        raise AssertionError("expected direct channel flux shape validation")


def test_aether_angular_reflective_boundary_supports_all_outer_faces():
    cases = [
        (
            (0, slice(None), slice(None)),
            (0, slice(None), slice(None), 0),
            (1, slice(None), slice(None), 0),
        ),
        (
            (-1, slice(None), slice(None)),
            (-1, slice(None), slice(None), 0),
            (-2, slice(None), slice(None), 0),
        ),
        (
            (slice(None), 0, slice(None)),
            (slice(None), 0, slice(None), 0),
            (slice(None), 1, slice(None), 0),
        ),
        (
            (slice(None), -1, slice(None)),
            (slice(None), -1, slice(None), 0),
            (slice(None), -2, slice(None), 0),
        ),
        (
            (slice(None), slice(None), 0),
            (slice(None), slice(None), 0, 0),
            (slice(None), slice(None), 1, 0),
        ),
        (
            (slice(None), slice(None), -1),
            (slice(None), slice(None), -1, 0),
            (slice(None), slice(None), -2, 0),
        ),
    ]

    for registration_key, boundary_slice, mirror_slice in cases:
        grid = fdtd.AetherGrid(shape=(4, 4, 4))
        response_rate_t = 20.0
        response_rate_p = 10.0
        grid.angular_momentum_t[:, :, :, 0] = 0.75
        grid.angular_momentum_p[:, :, :, 0] = 1.25
        grid.angular_momentum_t[boundary_slice] = 0.5
        grid.angular_momentum_p[boundary_slice] = 0.25
        grid[registration_key] = fdtd.AetherAngularReflectiveBoundary(
            response_rate_t=response_rate_t,
            response_rate_p=response_rate_p,
            sign_t=-1.0,
            sign_p=1.0,
        )

        source_t, source_p = grid.collect_native_angular_source_terms()

        expected_t = np.zeros((4, 4, 4, 1))
        expected_p = np.zeros((4, 4, 4, 1))
        expected_t[boundary_slice] = response_rate_t * (
            -grid.angular_momentum_t[mirror_slice]
            - grid.angular_momentum_t[boundary_slice]
        )
        expected_p[boundary_slice] = response_rate_p * (
            grid.angular_momentum_p[mirror_slice]
            - grid.angular_momentum_p[boundary_slice]
        )
        np.testing.assert_allclose(np.asarray(source_t), expected_t)
        np.testing.assert_allclose(np.asarray(source_p), expected_p)
        assert not np.any(np.asarray(grid.angular_tau))
        assert not np.any(np.asarray(grid.linear_a))


def test_aether_angular_reflective_boundary_rejects_non_outer_face():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid[1, :, :] = fdtd.AetherAngularReflectiveBoundary(
        response_rate_t=20.0,
        response_rate_p=10.0,
    )

    try:
        grid.collect_native_angular_source_terms()
    except ValueError as exc:
        assert "single outer grid face" in str(exc)
    else:
        raise AssertionError("expected outer-face validation")

    try:
        grid.boundaries[0].native_angular_boundary_contract()
    except ValueError as exc:
        assert "single outer grid face" in str(exc)
    else:
        raise AssertionError("expected outer-face contract validation")

    try:
        grid.collect_native_angular_boundary_contracts()
    except ValueError as exc:
        assert "single outer grid face" in str(exc)
    else:
        raise AssertionError("expected outer-face grid contract validation")

    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid[:, :, :] = fdtd.AetherAngularReflectiveBoundary(
        response_rate_t=20.0,
        response_rate_p=10.0,
    )

    try:
        grid.collect_native_angular_source_terms()
    except ValueError as exc:
        assert "single outer grid face" in str(exc)
    else:
        raise AssertionError("expected outer-face validation")

    try:
        grid.boundaries[0].native_angular_boundary_contract()
    except ValueError as exc:
        assert "single outer grid face" in str(exc)
    else:
        raise AssertionError("expected outer-face contract validation")


def test_aethergrid_passive_momentum_predictor_reflective_boundary_hook_64_step():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.05
    steps = 64
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.angular_momentum_t[0, :, :, 0] = 0.5
    grid.angular_momentum_p[0, :, :, 0] = 0.25
    initial_interior_t = np.asarray(grid.angular_momentum_t[1:]).copy()
    initial_interior_p = np.asarray(grid.angular_momentum_p[1:]).copy()
    expected_boundary_t = -np.asarray(grid.angular_momentum_t[1]).copy()
    expected_boundary_p = np.asarray(grid.angular_momentum_p[1]).copy()
    grid[0, :, :] = fdtd.AetherAngularReflectiveBoundary(
        response_rate_t=1.0 / delta,
        response_rate_p=1.0 / delta,
        sign_t=-1.0,
        sign_p=1.0,
    )

    for _ in range(steps):
        source_t, source_p = grid.collect_native_angular_source_terms()
        candidate_t, candidate_p = grid.predict_native_angular_momentum_step(
            delta=delta,
            source_t=source_t,
            source_p=source_p,
        )
        candidate_t = np.asarray(candidate_t)
        candidate_p = np.asarray(candidate_p)

        np.testing.assert_allclose(candidate_t[0], expected_boundary_t)
        np.testing.assert_allclose(candidate_p[0], expected_boundary_p)
        np.testing.assert_allclose(candidate_t[1:], initial_interior_t)
        np.testing.assert_allclose(candidate_p[1:], initial_interior_p)
        assert np.all(np.isfinite(candidate_t))
        assert np.all(np.isfinite(candidate_p))

        # Candidate-only reflection accounting benchmark; production dynamics
        # still do not call or apply the predictor.
        grid.angular_momentum_t = grid.native_angular_momentum_candidate_t
        grid.angular_momentum_p = grid.native_angular_momentum_candidate_p

    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_reflective_boundary_conserves_quadratic_energy_for_sign_flip():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.05
    steps = 64
    grid.angular_inertia_t[:, :, :, 0] = 1.0
    grid.angular_inertia_p[:, :, :, 0] = 1.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    _, _, initial_energy = grid.evaluate_native_angular_kinetic_energy()
    grid[0, :, :] = fdtd.AetherAngularReflectiveBoundary(
        response_rate_t=1.0 / delta,
        response_rate_p=1.0 / delta,
        sign_t=-1.0,
        sign_p=-1.0,
    )

    for _ in range(steps):
        source_t, source_p = grid.collect_native_angular_source_terms()
        candidate_t, candidate_p = grid.predict_native_angular_momentum_step(
            delta=delta,
            source_t=source_t,
            source_p=source_p,
        )
        _, _, candidate_energy = grid.evaluate_native_angular_kinetic_energy(
            momentum_t=candidate_t,
            momentum_p=candidate_p,
        )

        np.testing.assert_allclose(float(candidate_energy), float(initial_energy))
        assert np.all(np.isfinite(np.asarray(candidate_t)))
        assert np.all(np.isfinite(np.asarray(candidate_p)))

        # Candidate-only sign-reflection accounting benchmark; production
        # dynamics still do not call or apply the predictor.
        grid.angular_momentum_t = grid.native_angular_momentum_candidate_t
        grid.angular_momentum_p = grid.native_angular_momentum_candidate_p

    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_reflective_boundary_can_change_energy_when_unmatched():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.05
    grid.angular_inertia_t[:, :, :, 0] = 1.0
    grid.angular_inertia_p[:, :, :, 0] = 1.0
    grid.angular_momentum_t[:, :, :, 0] = 0.0
    grid.angular_momentum_p[:, :, :, 0] = 0.0
    grid.angular_momentum_t[0, :, :, 0] = 0.1
    grid.angular_momentum_p[0, :, :, 0] = 0.2
    grid.angular_momentum_t[1, :, :, 0] = 1.0
    grid.angular_momentum_p[1, :, :, 0] = 1.5
    _, _, initial_energy = grid.evaluate_native_angular_kinetic_energy()
    grid[0, :, :] = fdtd.AetherAngularReflectiveBoundary(
        response_rate_t=1.0 / delta,
        response_rate_p=1.0 / delta,
        sign_t=-1.0,
        sign_p=-1.0,
    )

    source_t, source_p = grid.collect_native_angular_source_terms()
    candidate_t, candidate_p = grid.predict_native_angular_momentum_step(
        delta=delta,
        source_t=source_t,
        source_p=source_p,
    )
    _, _, candidate_energy = grid.evaluate_native_angular_kinetic_energy(
        momentum_t=candidate_t,
        momentum_p=candidate_p,
    )

    assert float(candidate_energy) > float(initial_energy)
    assert np.all(np.isfinite(np.asarray(candidate_t)))
    assert np.all(np.isfinite(np.asarray(candidate_p)))
    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_transport_advance_reflective_boundary_conserves_quadratic_energy():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.05
    steps = 64
    grid.angular_inertia_t[:, :, :, 0] = 1.0
    grid.angular_inertia_p[:, :, :, 0] = 1.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    _, _, initial_energy = grid.evaluate_native_angular_kinetic_energy()
    grid[0, :, :] = fdtd.AetherAngularReflectiveBoundary(
        response_rate_t=1.0 / delta,
        response_rate_p=1.0 / delta,
        sign_t=-1.0,
        sign_p=-1.0,
    )

    for _ in range(steps):
        momentum_t, momentum_p = grid.advance_native_angular_momentum_transport(
            delta=delta,
            collect_sources=True,
        )
        _, _, current_energy = grid.evaluate_native_angular_kinetic_energy()

        np.testing.assert_allclose(float(current_energy), float(initial_energy))
        assert np.all(np.isfinite(np.asarray(momentum_t)))
        assert np.all(np.isfinite(np.asarray(momentum_p)))

    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.angular_A))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_step_reflective_boundary_conserves_quadratic_energy_64_step():
    grid = fdtd.AetherGrid(
        shape=(4, 4, 4),
        native_angular_transport=True,
    )
    steps = 64
    grid.angular_inertia_t[:, :, :, 0] = 1.0
    grid.angular_inertia_p[:, :, :, 0] = 1.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    _, _, initial_energy = grid.evaluate_native_angular_kinetic_energy()
    grid[0, :, :] = fdtd.AetherAngularReflectiveBoundary(
        response_rate_t=1.0 / grid.time_step,
        response_rate_p=1.0 / grid.time_step,
        sign_t=-1.0,
        sign_p=-1.0,
    )

    for _ in range(steps):
        grid.step()
        _, _, current_energy = grid.evaluate_native_angular_kinetic_energy()

        np.testing.assert_allclose(float(current_energy), float(initial_energy))
        assert np.all(np.isfinite(np.asarray(grid.angular_momentum_t)))
        assert np.all(np.isfinite(np.asarray(grid.angular_momentum_p)))
        assert np.all(np.isfinite(np.asarray(grid.linear_v)))
        assert np.all(np.isfinite(np.asarray(grid.angular_A)))

    assert grid.time_steps_passed == steps
    assert not np.any(np.asarray(grid.linear_v))


def test_aethergrid_passive_momentum_predictor_sponge_boundary_hook_64_step():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.05
    steps = 64
    damping_t = np.zeros((4, 4, 4, 1))
    damping_p = np.zeros((4, 4, 4, 1))
    damping_t[0, :, :, 0] = 0.25
    damping_t[-1, :, :, 0] = 0.25
    damping_p[0, :, :, 0] = 0.15
    damping_p[-1, :, :, 0] = 0.15
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
    metric_t = np.asarray(metric_t).copy()
    metric_p = np.asarray(metric_p).copy()
    grid[1, 1, 1] = NativeAngularExchangeProbe(metric_t, metric_p, "sources")
    grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )
    grid[-1, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()

    for step in range(steps):
        source_t, source_p = grid.collect_native_angular_source_terms()
        candidate_t, candidate_p = grid.predict_native_angular_momentum_step(
            delta=delta,
            source_t=source_t,
            source_p=source_p,
        )
        candidate_t = np.asarray(candidate_t)
        candidate_p = np.asarray(candidate_p)
        expected_t = initial_momentum_t * (1.0 - damping_t * delta) ** (step + 1)
        expected_p = initial_momentum_p * (1.0 - damping_p * delta) ** (step + 1)

        np.testing.assert_allclose(candidate_t, expected_t)
        np.testing.assert_allclose(candidate_p, expected_p)
        np.testing.assert_allclose(candidate_t[1:3], initial_momentum_t[1:3])
        np.testing.assert_allclose(candidate_p[1:3], initial_momentum_p[1:3])
        assert np.all(np.isfinite(candidate_t))
        assert np.all(np.isfinite(candidate_p))

        # Candidate-only boundary-hook accounting benchmark; production
        # dynamics still do not call or apply the predictor.
        grid.angular_momentum_t = grid.native_angular_momentum_candidate_t
        grid.angular_momentum_p = grid.native_angular_momentum_candidate_p

    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_sponge_boundary_reduces_quadratic_energy():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.05
    steps = 64
    damping_t = np.zeros((4, 4, 4, 1))
    damping_p = np.zeros((4, 4, 4, 1))
    damping_t[0, :, :, 0] = 0.25
    damping_t[-1, :, :, 0] = 0.25
    damping_p[0, :, :, 0] = 0.15
    damping_p[-1, :, :, 0] = 0.15
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()
    _, _, previous_energy = grid.evaluate_native_angular_kinetic_energy()
    grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )
    grid[-1, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )

    for step in range(steps):
        source_t, source_p = grid.collect_native_angular_source_terms()
        candidate_t, candidate_p = grid.predict_native_angular_momentum_step(
            delta=delta,
            source_t=source_t,
            source_p=source_p,
        )
        expected_t = initial_momentum_t * (1.0 - damping_t * delta) ** (step + 1)
        expected_p = initial_momentum_p * (1.0 - damping_p * delta) ** (step + 1)
        _, _, candidate_energy = grid.evaluate_native_angular_kinetic_energy(
            momentum_t=candidate_t,
            momentum_p=candidate_p,
        )
        expected_energy_t = 0.5 * expected_t * expected_t / 2.0
        expected_energy_p = 0.5 * expected_p * expected_p / 4.0
        expected_energy = np.sum(expected_energy_t + expected_energy_p)

        np.testing.assert_allclose(np.asarray(candidate_t), expected_t)
        np.testing.assert_allclose(np.asarray(candidate_p), expected_p)
        np.testing.assert_allclose(float(candidate_energy), expected_energy)
        assert float(candidate_energy) <= float(previous_energy)
        assert np.all(np.isfinite(np.asarray(candidate_t)))
        assert np.all(np.isfinite(np.asarray(candidate_p)))

        # Candidate-only sponge-energy accounting benchmark; production
        # dynamics still do not call or apply the predictor.
        grid.angular_momentum_t = grid.native_angular_momentum_candidate_t
        grid.angular_momentum_p = grid.native_angular_momentum_candidate_p
        previous_energy = candidate_energy

    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_transport_advance_sponge_boundary_reduces_quadratic_energy():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    delta = 0.05
    steps = 64
    damping_t = np.zeros((4, 4, 4, 1))
    damping_p = np.zeros((4, 4, 4, 1))
    damping_t[0, :, :, 0] = 0.25
    damping_t[-1, :, :, 0] = 0.25
    damping_p[0, :, :, 0] = 0.15
    damping_p[-1, :, :, 0] = 0.15
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()
    _, _, previous_energy = grid.evaluate_native_angular_kinetic_energy()
    grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )
    grid[-1, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )

    for step in range(steps):
        momentum_t, momentum_p = grid.advance_native_angular_momentum_transport(
            delta=delta,
            collect_sources=True,
        )
        expected_t = initial_momentum_t * (1.0 - damping_t * delta) ** (step + 1)
        expected_p = initial_momentum_p * (1.0 - damping_p * delta) ** (step + 1)
        _, _, current_energy = grid.evaluate_native_angular_kinetic_energy()
        expected_energy_t = 0.5 * expected_t * expected_t / 2.0
        expected_energy_p = 0.5 * expected_p * expected_p / 4.0
        expected_energy = np.sum(expected_energy_t + expected_energy_p)

        np.testing.assert_allclose(np.asarray(momentum_t), expected_t)
        np.testing.assert_allclose(np.asarray(momentum_p), expected_p)
        np.testing.assert_allclose(float(current_energy), expected_energy)
        assert float(current_energy) <= float(previous_energy)
        assert np.all(np.isfinite(np.asarray(momentum_t)))
        assert np.all(np.isfinite(np.asarray(momentum_p)))
        previous_energy = current_energy

    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.angular_A))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_step_sponge_boundary_reduces_quadratic_energy_64_step():
    grid = fdtd.AetherGrid(
        shape=(4, 4, 4),
        native_angular_transport=True,
    )
    steps = 64
    damping_t = np.zeros((4, 4, 4, 1))
    damping_p = np.zeros((4, 4, 4, 1))
    damping_t[0, :, :, 0] = 0.25
    damping_t[-1, :, :, 0] = 0.25
    damping_p[0, :, :, 0] = 0.15
    damping_p[-1, :, :, 0] = 0.15
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()
    _, _, previous_energy = grid.evaluate_native_angular_kinetic_energy()
    grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )
    grid[-1, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )

    for step in range(steps):
        grid.step()
        expected_t = (
            initial_momentum_t
            * (1.0 - damping_t * grid.time_step) ** (step + 1)
        )
        expected_p = (
            initial_momentum_p
            * (1.0 - damping_p * grid.time_step) ** (step + 1)
        )
        _, _, current_energy = grid.evaluate_native_angular_kinetic_energy()
        expected_energy_t = 0.5 * expected_t * expected_t / 2.0
        expected_energy_p = 0.5 * expected_p * expected_p / 4.0
        expected_energy = np.sum(expected_energy_t + expected_energy_p)

        np.testing.assert_allclose(np.asarray(grid.angular_momentum_t), expected_t)
        np.testing.assert_allclose(np.asarray(grid.angular_momentum_p), expected_p)
        np.testing.assert_allclose(float(current_energy), expected_energy)
        assert float(current_energy) <= float(previous_energy)
        assert np.all(np.isfinite(np.asarray(grid.angular_momentum_t)))
        assert np.all(np.isfinite(np.asarray(grid.angular_momentum_p)))
        assert np.all(np.isfinite(np.asarray(grid.linear_v)))
        previous_energy = current_energy

    assert grid.time_steps_passed == steps
    assert not np.any(np.asarray(grid.linear_v))


def test_aethergrid_step_sponge_boundary_matches_exchange_power_balance():
    grid = fdtd.AetherGrid(
        shape=(4, 4, 4),
        grid_spacing=1.0e6,
        native_angular_transport=True,
    )
    steps = 64
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )
    grid[-1, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )

    for _ in range(steps):
        _, _, energy_before = grid.evaluate_native_angular_kinetic_energy()
        source_t, source_p = grid.collect_native_angular_source_terms()
        source_t = np.asarray(source_t).copy()
        source_p = np.asarray(source_p).copy()
        _, _, exchange_power = grid.evaluate_native_angular_exchange_power(
            source_t=source_t,
            source_p=source_p,
        )
        discrete_correction = 0.5 * grid.time_step**2 * np.sum(
            source_t * source_t / np.asarray(grid.angular_inertia_t)
            + source_p * source_p / np.asarray(grid.angular_inertia_p)
        )

        grid.step()
        _, _, energy_after = grid.evaluate_native_angular_kinetic_energy()

        np.testing.assert_allclose(
            float(energy_after) - float(energy_before),
            grid.time_step * float(exchange_power) + discrete_correction,
        )
        assert float(exchange_power) <= 0.0
        assert float(energy_after) <= float(energy_before)
        assert np.all(np.isfinite(np.asarray(grid.native_angular_exchange_power_t)))
        assert np.all(np.isfinite(np.asarray(grid.native_angular_exchange_power_p)))

    assert grid.time_steps_passed == steps
    assert not np.any(np.asarray(grid.linear_v))


def test_aethergrid_step_dynamic_transport_sponge_boundary_is_dissipative():
    def run_dynamic_transport(use_sponge):
        grid = fdtd.AetherGrid(
            shape=(5, 5, 5),
            grid_spacing=8.25e6,
            native_angular_transport=True,
            native_angular_collect_sources=True,
            native_angular_update_clocks=True,
        )
        steps = 64
        grid.ell_t[:, :, :, 0] = 0.1
        grid.ell_p[:, :, :, 0] = 0.1
        grid.angular_inertia_t[:, :, :, 0] = 1.0
        grid.angular_inertia_p[:, :, :, 0] = 1.0
        initial_rate_t = np.zeros((5, 5, 5, 1))
        initial_rate_p = np.zeros((5, 5, 5, 1))
        initial_rate_t[2, 2, 2, 0] = 0.5
        initial_rate_p[2, 2, 2, 0] = -0.25
        previous_t = np.asarray(grid.angular_momentum_t).copy()
        previous_p = np.asarray(grid.angular_momentum_p).copy()
        previous_t -= grid.time_step * initial_rate_t
        previous_p -= grid.time_step * initial_rate_p
        energies = []
        boundary_momenta = []
        exchange_powers = []

        if use_sponge:
            grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
                damping_t=0.5,
                damping_p=0.5,
            )
            grid[-1, :, :] = fdtd.AetherAngularSpongeBoundary(
                damping_t=0.5,
                damping_p=0.5,
            )

        for _ in range(steps):
            current_t = np.asarray(grid.angular_momentum_t).copy()
            current_p = np.asarray(grid.angular_momentum_p).copy()
            grid.update_native_angular_torque(
                previous_momentum_t=previous_t,
                previous_momentum_p=previous_p,
                delta=grid.time_step,
            )
            grid.project_native_angular_torque()
            source_t, source_p = grid.collect_native_angular_source_terms()
            _, _, exchange_power = grid.evaluate_native_angular_exchange_power(
                source_t=source_t,
                source_p=source_p,
            )

            grid.step()

            _, _, energy = grid.evaluate_native_angular_kinetic_energy()
            momentum_t = np.asarray(grid.angular_momentum_t)
            momentum_p = np.asarray(grid.angular_momentum_p)
            energies.append(float(energy))
            boundary_momenta.append(
                float(np.sum(momentum_t[[0, -1]] ** 2 + momentum_p[[0, -1]] ** 2))
            )
            exchange_powers.append(float(exchange_power))
            assert np.all(np.isfinite(momentum_t))
            assert np.all(np.isfinite(momentum_p))

            previous_t = current_t
            previous_p = current_p

        assert grid.time_steps_passed == steps
        assert not np.any(np.asarray(grid.linear_v))
        return (
            np.asarray(energies),
            np.asarray(boundary_momenta),
            np.asarray(exchange_powers),
        )

    closed_energy, closed_boundary, closed_power = run_dynamic_transport(False)
    sponge_energy, sponge_boundary, sponge_power = run_dynamic_transport(True)

    np.testing.assert_allclose(closed_power, 0.0)
    assert np.all(sponge_power <= 0.0)
    assert np.any(sponge_power < 0.0)
    assert np.all(sponge_energy <= closed_energy)
    assert np.all(sponge_boundary <= closed_boundary)
    assert sponge_energy[-1] < closed_energy[-1]
    assert sponge_boundary[-1] < closed_boundary[-1]


def test_aethergrid_step_static_torque_transport_matches_power_balance():
    grid = fdtd.AetherGrid(
        shape=(4, 4, 4),
        grid_spacing=1.0e6,
        native_angular_transport=True,
        native_angular_collect_sources=False,
    )
    steps = 64
    grid.ell_t[:, :, :, 0] = 0.5
    grid.ell_p[:, :, :, 0] = 0.75
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.native_angular_tau[1, 1, 1, 0] = 0.2
    grid.native_angular_tau[2, 1, 1, 1] = -0.1
    _, _, initial_energy = grid.evaluate_native_angular_kinetic_energy()
    max_energy = float(initial_energy)

    for _ in range(steps):
        _, _, energy_before = grid.evaluate_native_angular_kinetic_energy()
        _, _, transport_power = grid.evaluate_native_angular_transport_power()
        rhs_t = np.asarray(grid.native_angular_momentum_rhs_t).copy()
        rhs_p = np.asarray(grid.native_angular_momentum_rhs_p).copy()
        discrete_correction = 0.5 * grid.time_step**2 * np.sum(
            rhs_t * rhs_t / np.asarray(grid.angular_inertia_t)
            + rhs_p * rhs_p / np.asarray(grid.angular_inertia_p)
        )

        grid.step()
        _, _, energy_after = grid.evaluate_native_angular_kinetic_energy()

        np.testing.assert_allclose(
            float(energy_after) - float(energy_before),
            grid.time_step * float(transport_power) + discrete_correction,
        )
        assert np.all(np.isfinite(np.asarray(grid.angular_momentum_t)))
        assert np.all(np.isfinite(np.asarray(grid.angular_momentum_p)))
        assert np.all(np.isfinite(np.asarray(grid.native_angular_transport_power_t)))
        assert np.all(np.isfinite(np.asarray(grid.native_angular_transport_power_p)))
        max_energy = max(max_energy, float(energy_after))

    assert max_energy < 1.01 * float(initial_energy)
    assert grid.time_steps_passed == steps
    assert not np.any(np.asarray(grid.native_angular_source_t))
    assert not np.any(np.asarray(grid.native_angular_source_p))
    assert not np.any(np.asarray(grid.linear_v))


def test_aethergrid_step_native_angular_transport_updates_clocks_64_step():
    grid = fdtd.AetherGrid(
        shape=(4, 4, 4),
        grid_spacing=1.0e6,
        native_angular_transport=True,
        native_angular_update_clocks=True,
    )
    steps = 64
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )
    grid[-1, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )

    for _ in range(steps):
        grid.step()
        np.testing.assert_allclose(
            np.asarray(grid.omega_t),
            np.asarray(grid.angular_momentum_t) / np.asarray(grid.angular_inertia_t),
        )
        np.testing.assert_allclose(
            np.asarray(grid.omega_p),
            np.asarray(grid.angular_momentum_p) / np.asarray(grid.angular_inertia_p),
        )
        assert np.all(np.isfinite(np.asarray(grid.omega_t)))
        assert np.all(np.isfinite(np.asarray(grid.omega_p)))
        assert np.all(np.isfinite(np.asarray(grid.angular_momentum_t)))
        assert np.all(np.isfinite(np.asarray(grid.angular_momentum_p)))

    assert grid.time_steps_passed == steps
    assert not np.any(np.asarray(grid.linear_v))


def test_aethergrid_step_native_angular_clock_update_requires_positive_inertia():
    grid = fdtd.AetherGrid(
        shape=(2, 2, 2),
        native_angular_transport=True,
        native_angular_update_clocks=True,
    )
    grid.angular_inertia_t[:, :, :, 0] = 1.0
    grid.angular_inertia_p[:, :, :, 0] = 0.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()

    try:
        grid.step()
    except ValueError as exc:
        assert "inertia" in str(exc)
    else:
        raise AssertionError("expected positive-inertia validation")

    assert grid.time_steps_passed == 0
    np.testing.assert_allclose(np.asarray(grid.angular_momentum_t), initial_momentum_t)
    np.testing.assert_allclose(np.asarray(grid.angular_momentum_p), initial_momentum_p)


def test_aethergrid_transport_residual_matches_staged_torque_minus_rhs():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_torque_t[:, :, :, 0] = 0.25
    grid.angular_torque_p[:, :, :, 0] = 0.5
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    source_t = np.ones((4, 4, 4, 1)) * 0.1
    source_p = np.ones((4, 4, 4, 1)) * 0.2

    residual_t, residual_p = grid.evaluate_native_angular_transport_residual(
        source_t=source_t,
        source_p=source_p,
    )

    np.testing.assert_allclose(
        np.asarray(residual_t),
        np.asarray(grid.angular_torque_t) - np.asarray(grid.native_angular_momentum_rhs_t),
    )
    np.testing.assert_allclose(
        np.asarray(residual_p),
        np.asarray(grid.angular_torque_p) - np.asarray(grid.native_angular_momentum_rhs_p),
    )


def test_aethergrid_transport_residual_is_bounded_and_passive_when_repeated():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.ell_t[:, :, :, 0] = 1.5
    grid.ell_p[:, :, :, 0] = 2.5
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid.angular_torque_t[:, :, :, 0] = 0.25
    grid.angular_torque_p[:, :, :, 0] = 0.5
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    initial_momentum_t = np.asarray(grid.angular_momentum_t).copy()
    initial_momentum_p = np.asarray(grid.angular_momentum_p).copy()

    residual_t, residual_p = grid.evaluate_native_angular_transport_residual()
    first_t = np.asarray(residual_t).copy()
    first_p = np.asarray(residual_p).copy()

    for _ in range(20):
        residual_t, residual_p = grid.evaluate_native_angular_transport_residual()
        np.testing.assert_allclose(np.asarray(residual_t), first_t)
        np.testing.assert_allclose(np.asarray(residual_p), first_p)
        assert np.all(np.isfinite(np.asarray(residual_t)))
        assert np.all(np.isfinite(np.asarray(residual_p)))
        assert np.all(np.isfinite(np.asarray(grid.native_angular_momentum_rhs_t)))
        assert np.all(np.isfinite(np.asarray(grid.native_angular_momentum_rhs_p)))

    assert grid.time_steps_passed == 0
    np.testing.assert_allclose(np.asarray(grid.angular_momentum_t), initial_momentum_t)
    np.testing.assert_allclose(np.asarray(grid.angular_momentum_p), initial_momentum_p)
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_collects_explicit_native_angular_source_terms():
    grid = fdtd.AetherGrid(shape=(3, 3, 3))
    source_t = np.ones((3, 3, 3, 1)) * 0.25
    source_p = np.ones((3, 3, 3, 1)) * 0.5
    boundary_t = np.ones((3, 3, 3, 1)) * 0.75
    boundary_p = np.ones((3, 3, 3, 1)) * 1.25

    grid[1, 1, 1] = NativeAngularExchangeProbe(source_t, source_p, "sources")
    grid[0, :, :] = NativeAngularExchangeProbe(
        boundary_t,
        boundary_p,
        "boundaries",
    )

    collected_t, collected_p = grid.collect_native_angular_source_terms()

    np.testing.assert_allclose(np.asarray(collected_t), 1.0)
    np.testing.assert_allclose(np.asarray(collected_p), 1.75)
    assert collected_t is grid.native_angular_source_t
    assert collected_p is grid.native_angular_source_p
    assert grid.time_steps_passed == 0
    assert not np.any(np.asarray(grid.angular_tau))
    assert not np.any(np.asarray(grid.linear_a))


def test_aethergrid_collected_native_angular_sources_can_drive_residual():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.ell_t[:, :, :, 0] = 2.0
    grid.ell_p[:, :, :, 0] = 3.0
    grid.angular_torque_t[:, :, :, 0] = 0.25
    grid.angular_torque_p[:, :, :, 0] = 0.5
    grid.native_angular_tau[1, 1, 1, 0] = 2.0
    grid.native_angular_tau[2, 1, 1, 1] = -1.0
    metric_t, metric_p = grid.evaluate_native_angular_metric_divergence()
    source_t = np.asarray(grid.angular_torque_t) + np.asarray(metric_t)
    source_p = np.asarray(grid.angular_torque_p) + np.asarray(metric_p)
    grid[0, :, :] = NativeAngularExchangeProbe(
        source_t,
        source_p,
        "boundaries",
    )

    collected_t, collected_p = grid.collect_native_angular_source_terms()
    residual_t, residual_p = grid.evaluate_native_angular_transport_residual(
        source_t=collected_t,
        source_p=collected_p,
    )

    np.testing.assert_allclose(np.asarray(residual_t), 0.0)
    np.testing.assert_allclose(np.asarray(residual_p), 0.0)


def test_aethergrid_collects_only_requested_native_angular_sources():
    grid = fdtd.AetherGrid(shape=(3, 3, 3))
    source_t = np.ones((3, 3, 3, 1)) * 0.25
    source_p = np.ones((3, 3, 3, 1)) * 0.5
    boundary_t = np.ones((3, 3, 3, 1)) * 0.75
    boundary_p = np.ones((3, 3, 3, 1)) * 1.25

    grid[1, 1, 1] = NativeAngularExchangeProbe(source_t, source_p, "sources")
    grid[0, :, :] = NativeAngularExchangeProbe(
        boundary_t,
        boundary_p,
        "boundaries",
    )

    collected_t, collected_p = grid.collect_native_angular_source_terms(
        include_boundaries=False,
    )

    np.testing.assert_allclose(np.asarray(collected_t), 0.25)
    np.testing.assert_allclose(np.asarray(collected_p), 0.5)


def test_aethergrid_native_angular_source_collection_resets_each_call():
    grid = fdtd.AetherGrid(shape=(3, 3, 3))
    source_t = np.ones((3, 3, 3, 1)) * 0.25
    source_p = np.ones((3, 3, 3, 1)) * 0.5
    boundary_t = np.ones((3, 3, 3, 1)) * 0.75
    boundary_p = np.ones((3, 3, 3, 1)) * 1.25

    grid[1, 1, 1] = NativeAngularExchangeProbe(source_t, source_p, "sources")
    grid[0, :, :] = NativeAngularExchangeProbe(
        boundary_t,
        boundary_p,
        "boundaries",
    )

    collected_t, collected_p = grid.collect_native_angular_source_terms()
    np.testing.assert_allclose(np.asarray(collected_t), 1.0)
    np.testing.assert_allclose(np.asarray(collected_p), 1.75)

    grid.native_angular_source_t[:, :, :, 0] = 99.0
    grid.native_angular_source_p[:, :, :, 0] = 101.0
    collected_t, collected_p = grid.collect_native_angular_source_terms(
        include_sources=False,
        include_boundaries=False,
    )
    np.testing.assert_allclose(np.asarray(collected_t), 0.0)
    np.testing.assert_allclose(np.asarray(collected_p), 0.0)

    collected_t, collected_p = grid.collect_native_angular_source_terms(
        include_sources=False,
    )
    np.testing.assert_allclose(np.asarray(collected_t), 0.75)
    np.testing.assert_allclose(np.asarray(collected_p), 1.25)

    collected_t, collected_p = grid.collect_native_angular_source_terms()
    np.testing.assert_allclose(np.asarray(collected_t), 1.0)
    np.testing.assert_allclose(np.asarray(collected_p), 1.75)
    collected_t, collected_p = grid.collect_native_angular_source_terms()
    np.testing.assert_allclose(np.asarray(collected_t), 1.0)
    np.testing.assert_allclose(np.asarray(collected_p), 1.75)


def test_aethergrid_rejects_malformed_native_angular_source_terms():
    grid = fdtd.AetherGrid(shape=(3, 3, 3))
    source_t = np.ones((3, 3, 3, 1))
    source_p = np.ones((3, 3, 3))
    grid[1, 1, 1] = NativeAngularExchangeProbe(source_t, source_p, "sources")

    try:
        grid.collect_native_angular_source_terms()
    except ValueError as exc:
        assert "native angular poloidal source" in str(exc)
    else:
        raise AssertionError("expected ValueError for malformed source shape")


def test_aethergrid_step_does_not_promote_cartesian_omega_to_native_clocks():
    grid = fdtd.AetherGrid(shape=(5, 5, 5))
    grid[2, 2, 2] = fdtd.AetherPointSource(amplitude=1.0, phase_shift=pi / 2)

    grid.step()

    assert np.any(np.asarray(grid.angular_omega) != 0.0)
    assert not np.any(np.asarray(grid.omega_t))
    assert not np.any(np.asarray(grid.omega_p))
    assert not np.any(np.asarray(grid.angular_gamma))
    assert not np.any(np.asarray(grid.angular_inertia_t))
    assert not np.any(np.asarray(grid.angular_inertia_p))
    assert not np.any(np.asarray(grid.angular_momentum_t))
    assert not np.any(np.asarray(grid.angular_momentum_p))
    assert not np.any(np.asarray(grid.angular_torque_t))
    assert not np.any(np.asarray(grid.angular_torque_p))
    assert not np.any(np.asarray(grid.native_angular_tau))
    assert not np.any(np.asarray(grid.native_angular_tau_divergence))
    assert not np.any(np.asarray(grid.native_angular_tau_metric_divergence_t))
    assert not np.any(np.asarray(grid.native_angular_tau_metric_divergence_p))
    assert not np.any(np.asarray(grid.native_angular_transport_residual_t))
    assert not np.any(np.asarray(grid.native_angular_transport_residual_p))
    assert not np.any(np.asarray(grid.native_angular_source_t))
    assert not np.any(np.asarray(grid.native_angular_source_p))
    assert not np.any(np.asarray(grid.native_angular_momentum_rhs_t))
    assert not np.any(np.asarray(grid.native_angular_momentum_rhs_p))
    assert not np.any(np.asarray(grid.native_angular_momentum_candidate_t))
    assert not np.any(np.asarray(grid.native_angular_momentum_candidate_p))


def test_aethergrid_reset_preserves_native_angular_geometry():
    grid = fdtd.AetherGrid(shape=(3, 3, 3))
    grid.ell_t[1, 1, 1, 0] = 2.0
    grid.ell_p[1, 1, 1, 0] = 3.0
    grid.angular_inertia_t[:, :, :, 0] = 4.0
    grid.angular_inertia_p[:, :, :, 0] = 5.0
    grid.omega_t[1, 1, 1, 0] = 1.0
    grid.evaluate_angular_clock_benchmark(delta=0.2)
    grid.update_native_angular_momentum()
    grid.update_native_angular_torque(delta=0.2)
    grid.project_native_angular_torque()
    grid.evaluate_native_angular_torque_divergence()
    grid.evaluate_native_angular_metric_divergence()
    grid.evaluate_native_angular_transport_residual()
    grid.predict_native_angular_momentum_step(delta=0.2)
    grid.native_angular_source_t[1, 1, 1, 0] = 6.0
    grid.native_angular_source_p[1, 1, 1, 0] = 7.0
    grid.evaluate_native_angular_exchange_power()
    grid.evaluate_native_angular_transport_power()
    grid.evaluate_native_angular_linear_response()
    grid.linear_v[1, 1, 1, 0] = 1.0
    grid.evaluate_linear_charge_flux_candidate(delta_t=0.2)
    grid.evaluate_native_angular_charge_flux_candidates(delta_t=0.2)
    grid.reduce_native_angular_charge_flux(mode="additive")
    grid.evaluate_native_angular_boundary_channel_flux(axis="x", side="low")
    grid.evaluate_native_angular_boundary_direct_channel_flux(axis="x", side="low")

    grid.reset()

    assert float(grid.ell_t[1, 1, 1, 0]) == 2.0
    assert float(grid.ell_p[1, 1, 1, 0]) == 3.0
    assert float(grid.angular_inertia_t[1, 1, 1, 0]) == 4.0
    assert float(grid.angular_inertia_p[1, 1, 1, 0]) == 5.0
    assert not np.any(np.asarray(grid.omega_t))
    assert not np.any(np.asarray(grid.angular_clock_lambda))
    assert not np.any(np.asarray(grid.angular_momentum_t))
    assert not np.any(np.asarray(grid.angular_momentum_p))
    assert not np.any(np.asarray(grid.angular_torque_t))
    assert not np.any(np.asarray(grid.angular_torque_p))
    assert not np.any(np.asarray(grid.native_angular_tau))
    assert not np.any(np.asarray(grid.native_angular_tau_divergence))
    assert not np.any(np.asarray(grid.native_angular_tau_metric_divergence_t))
    assert not np.any(np.asarray(grid.native_angular_tau_metric_divergence_p))
    assert not np.any(np.asarray(grid.native_angular_transport_residual_t))
    assert not np.any(np.asarray(grid.native_angular_transport_residual_p))
    assert not np.any(np.asarray(grid.native_angular_source_t))
    assert not np.any(np.asarray(grid.native_angular_source_p))
    assert not np.any(np.asarray(grid.native_angular_exchange_power_t))
    assert not np.any(np.asarray(grid.native_angular_exchange_power_p))
    assert not np.any(np.asarray(grid.native_angular_momentum_rhs_t))
    assert not np.any(np.asarray(grid.native_angular_momentum_rhs_p))
    assert not np.any(np.asarray(grid.native_angular_transport_power_t))
    assert not np.any(np.asarray(grid.native_angular_transport_power_p))
    assert not np.any(np.asarray(grid.native_angular_momentum_candidate_t))
    assert not np.any(np.asarray(grid.native_angular_momentum_candidate_p))
    assert not np.any(np.asarray(grid.native_angular_linear_response))
    assert not np.any(np.asarray(grid.linear_charge_flux_candidate))
    assert not np.any(np.asarray(grid.native_angular_charge_flux_t))
    assert not np.any(np.asarray(grid.native_angular_charge_flux_p))
    assert not np.any(np.asarray(grid.native_angular_charge_reduction))
    assert not np.any(np.asarray(grid.native_angular_boundary_normal_flux))
    assert not np.any(np.asarray(grid.native_angular_boundary_incident_flux))
    assert not np.any(np.asarray(grid.native_angular_boundary_outgoing_flux))
    assert not np.any(np.asarray(grid.native_angular_boundary_incident_flux_t))
    assert not np.any(np.asarray(grid.native_angular_boundary_incident_flux_p))
    assert not np.any(np.asarray(grid.native_angular_boundary_outgoing_flux_t))
    assert not np.any(np.asarray(grid.native_angular_boundary_outgoing_flux_p))
    assert not np.any(
        np.asarray(grid.native_angular_boundary_direct_normal_flux_t)
    )
    assert not np.any(
        np.asarray(grid.native_angular_boundary_direct_normal_flux_p)
    )
    assert not np.any(
        np.asarray(grid.native_angular_boundary_direct_incident_flux_t)
    )
    assert not np.any(
        np.asarray(grid.native_angular_boundary_direct_incident_flux_p)
    )
    assert not np.any(
        np.asarray(grid.native_angular_boundary_direct_outgoing_flux_t)
    )
    assert not np.any(
        np.asarray(grid.native_angular_boundary_direct_outgoing_flux_p)
    )


def test_aether_point_source_injects_native_velocity():
    grid = fdtd.AetherGrid(shape=(5, 5, 5))
    grid[2, 2, 2] = fdtd.AetherPointSource(amplitude=1.0, phase_shift=pi / 2)

    grid.step()

    assert float(grid.v[2, 2, 2, 2]) > 0.0
    assert np.any(np.asarray(grid.omega) != 0.0)


def test_aether_line_source_injects_native_velocity_profile():
    grid = fdtd.AetherGrid(shape=(6, 6, 6))
    grid[1:5, 1:5, 1:5] = fdtd.AetherLineSource(
        amplitude=1.0,
        phase_shift=pi / 2,
    )

    grid.step()

    assert np.any(np.asarray(grid.v[..., 2]) > 0.0)


def test_aether_source_accepts_explicit_linear_field_name():
    grid = fdtd.AetherGrid(shape=(5, 5, 5))
    grid[2, 2, 2] = fdtd.AetherPointSource(
        amplitude=1.0,
        phase_shift=pi / 2,
        field="linear_v",
    )

    grid.step()

    assert float(grid.linear_v[2, 2, 2, 2]) > 0.0
    assert float(grid.v[2, 2, 2, 2]) > 0.0


def test_aether_native_angular_point_source_exposes_source_terms():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid[1, 1, 1] = fdtd.AetherNativeAngularPointSource(
        amplitude_t=0.5,
        amplitude_p=1.25,
        phase_shift=pi / 2,
        name="native_source",
    )

    source_t, source_p = grid.collect_native_angular_source_terms()

    expected_t = np.zeros((4, 4, 4, 1))
    expected_p = np.zeros((4, 4, 4, 1))
    expected_t[1, 1, 1, 0] = 0.5
    expected_p[1, 1, 1, 0] = 1.25
    np.testing.assert_allclose(np.asarray(source_t), expected_t)
    np.testing.assert_allclose(np.asarray(source_p), expected_p)
    assert grid.native_source is grid.sources[0]
    assert not np.any(np.asarray(grid.angular_momentum_t))
    assert not np.any(np.asarray(grid.angular_momentum_p))
    assert grid.time_steps_passed == 0


def test_aether_native_angular_point_source_drives_detected_transport_step():
    grid = fdtd.AetherGrid(
        shape=(4, 4, 4),
        grid_spacing=1.0e6,
        native_angular_transport=True,
        native_angular_update_clocks=True,
    )
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid[1, 1, 1] = fdtd.AetherNativeAngularPointSource(
        amplitude_t=0.5,
        amplitude_p=1.25,
        phase_shift=pi / 2,
    )
    grid[1:1, 1:1, 1:1] = fdtd.AetherNativeAngularDetector(
        name="native_detector",
        fields=(
            "angular_momentum_t",
            "angular_momentum_p",
            "omega_t",
            "omega_p",
            "native_angular_source_t",
            "native_angular_source_p",
        ),
        record_energy=True,
    )

    grid.step()

    expected_momentum_t = grid.time_step * 0.5
    expected_momentum_p = grid.time_step * 1.25
    np.testing.assert_allclose(
        float(grid.angular_momentum_t[1, 1, 1, 0]),
        expected_momentum_t,
    )
    np.testing.assert_allclose(
        float(grid.angular_momentum_p[1, 1, 1, 0]),
        expected_momentum_p,
    )
    detected_t = np.asarray(grid.native_detector.readings["angular_momentum_t"][0])
    detected_p = np.asarray(grid.native_detector.readings["angular_momentum_p"][0])
    detected_source_t = np.asarray(
        grid.native_detector.readings["native_angular_source_t"][0]
    )
    detected_source_p = np.asarray(
        grid.native_detector.readings["native_angular_source_p"][0]
    )
    np.testing.assert_allclose(detected_t[0, 0, 0, 0], expected_momentum_t)
    np.testing.assert_allclose(detected_p[0, 0, 0, 0], expected_momentum_p)
    np.testing.assert_allclose(detected_source_t[0, 0, 0, 0], 0.5)
    np.testing.assert_allclose(detected_source_p[0, 0, 0, 0], 1.25)
    np.testing.assert_allclose(float(grid.omega_t[1, 1, 1, 0]), expected_momentum_t / 2.0)
    np.testing.assert_allclose(float(grid.omega_p[1, 1, 1, 0]), expected_momentum_p / 4.0)
    assert len(grid.native_detector.energy) == 1
    assert float(grid.native_detector.energy[0]) > 0.0
    assert grid.time_steps_passed == 1
    assert not np.any(np.asarray(grid.linear_v))


def test_aether_native_angular_observable_example_runs():
    package_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = "."

    result = subprocess.run(
        [
            sys.executable,
            "examples/aether_native_angular_observables.py",
            "--backend",
            "numpy",
            "--steps",
            "8",
        ],
        cwd=package_root,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "Native-angular observable scene" in result.stdout
    assert "steps: 8" in result.stdout
    assert "detector_samples: 8" in result.stdout
    assert "finite_native_momentum: True" in result.stdout
    assert "finite_native_clocks: True" in result.stdout
    assert "linear_v_quiet: True" in result.stdout


def test_aether_native_angular_dynamic_transport_example_runs():
    package_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = "."

    result = subprocess.run(
        [
            sys.executable,
            "examples/aether_native_angular_dynamic_transport.py",
            "--backend",
            "numpy",
            "--steps",
            "16",
            "--sponge",
            "--feedback",
            "--charge-observables",
        ],
        cwd=package_root,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "Dynamic native-angular transport scene" in result.stdout
    assert "steps: 16" in result.stdout
    assert "sponge_enabled: True" in result.stdout
    assert "feedback_enabled: True" in result.stdout
    assert "detector_samples: 16" in result.stdout
    assert "finite_native_momentum: True" in result.stdout
    assert "finite_linear_v: True" in result.stdout
    assert "nonzero_transport_seen: True" in result.stdout
    assert "boundary_contracts: 2" in result.stdout
    assert "boundary_contracts_source_accounting: True" in result.stdout
    assert "physical_boundary_laws: False" in result.stdout
    assert "max_boundary_incident_flux:" in result.stdout
    assert "max_boundary_outgoing_flux:" in result.stdout
    assert "max_abs_boundary_net_flux:" in result.stdout
    assert "charge_observables_enabled: True" in result.stdout
    assert "finite_charge_candidates: True" in result.stdout


def test_aether_native_angular_boundary_flux_example_runs():
    package_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = "."

    result = subprocess.run(
        [
            sys.executable,
            "examples/aether_native_angular_boundary_flux.py",
            "--backend",
            "numpy",
            "--steps",
            "32",
        ],
        cwd=package_root,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "Native-angular boundary-flux acceptance scene" in result.stdout
    assert "steps: 32" in result.stdout
    assert "no_exchange_finite: True" in result.stdout
    assert "sponge_finite: True" in result.stdout
    assert "matched_flux_finite: True" in result.stdout
    assert "direct_matched_flux_finite: True" in result.stdout
    assert "matched_flux_swapped_frame_finite: True" in result.stdout
    assert "no_exchange_boundary_contracts: 2" in result.stdout
    assert "sponge_boundary_contracts: 2" in result.stdout
    assert "matched_flux_boundary_contracts: 2" in result.stdout
    assert "direct_matched_flux_boundary_contracts: 2" in result.stdout
    assert "matched_flux_swapped_frame_boundary_contracts: 2" in result.stdout
    assert "no_exchange_physical_boundary_laws: False" in result.stdout
    assert "sponge_physical_boundary_laws: False" in result.stdout
    assert "matched_flux_physical_boundary_laws: True" in result.stdout
    assert "direct_matched_flux_physical_boundary_laws: True" in result.stdout
    assert "matched_flux_swapped_frame_physical_boundary_laws: True" in result.stdout
    assert "matched_flux_max_incident_flux_t:" in result.stdout
    assert "matched_flux_max_incident_flux_p:" in result.stdout
    assert "matched_flux_swapped_frame_max_incident_flux_t:" in result.stdout
    assert "matched_flux_swapped_frame_max_incident_flux_p:" in result.stdout
    assert "matched_flux_xy_max_incident_flux_t:" in result.stdout
    assert "matched_flux_xy_max_incident_flux_p:" in result.stdout
    assert "matched_flux_max_outgoing_flux_t:" in result.stdout
    assert "matched_flux_max_outgoing_flux_p:" in result.stdout
    assert "matched_flux_swapped_frame_max_outgoing_flux_t:" in result.stdout
    assert "matched_flux_swapped_frame_max_outgoing_flux_p:" in result.stdout
    assert "matched_flux_xy_max_outgoing_flux_t:" in result.stdout
    assert "matched_flux_xy_max_outgoing_flux_p:" in result.stdout
    assert "sponge_energy_no_greater: True" in result.stdout
    assert "matched_flux_energy_no_greater: True" in result.stdout
    assert "direct_matched_flux_energy_no_greater: True" in result.stdout
    assert "matched_flux_swapped_frame_energy_no_greater: True" in result.stdout
    assert "matched_flux_xy_energy_no_greater: True" in result.stdout
    assert "sponge_boundary_momentum_no_greater: True" in result.stdout
    assert "matched_flux_boundary_momentum_no_greater: True" in result.stdout
    assert "direct_matched_flux_boundary_momentum_no_greater: True" in result.stdout
    assert "matched_flux_swapped_frame_boundary_momentum_no_greater:" in result.stdout
    assert "matched_flux_xy_boundary_momentum_no_greater: True" in result.stdout
    assert "sponge_exchange_nonpositive: True" in result.stdout
    assert "matched_flux_exchange_nonpositive: True" in result.stdout
    assert "direct_matched_flux_exchange_nonpositive: True" in result.stdout
    assert "matched_flux_swapped_frame_exchange_nonpositive: True" in result.stdout
    assert "matched_flux_xy_exchange_nonpositive: True" in result.stdout
    assert "counterpropagating_direct_exposes_hidden_flux: True" in result.stdout
    assert "counterpropagating_matched_no_exchange: True" in result.stdout
    assert (
        "counterpropagating_direct_matched_absorbs_hidden_flux: True"
        in result.stdout
    )
    assert (
        "counterpropagating_propagated_hidden_flux_at_launch: True"
        in result.stdout
    )
    assert (
        "counterpropagating_propagated_matched_no_launch_exchange: True"
        in result.stdout
    )
    assert (
        "counterpropagating_propagated_direct_absorbs_launch_hidden_flux: True"
        in result.stdout
    )
    assert (
        "counterpropagating_propagated_direct_energy_below_matched: True"
        in result.stdout
    )
    assert (
        "counterpropagating_propagated_direct_boundary_momentum_below_matched: True"
        in result.stdout
    )
    assert (
        "counterpropagating_propagated_direct_no_artificial_momentum_injection: True"
        in result.stdout
    )
    assert (
        "counterpropagating_propagated_decision: "
        "direct_channel_candidate_preferred_for_hidden_flux_benchmark"
        in result.stdout
    )
    summary = {
        key: value
        for key, value in (
            line.strip().split(": ", 1)
            for line in result.stdout.splitlines()
            if ": " in line
        )
    }
    assert float(summary["matched_flux_max_incident_flux_p"]) == 0.0
    assert float(summary["matched_flux_swapped_frame_max_incident_flux_t"]) == 0.0
    assert float(summary["matched_flux_swapped_frame_max_incident_flux_p"]) > 0.0
    assert float(summary["matched_flux_swapped_frame_final_energy"]) < float(
        summary["no_exchange_final_energy"]
    )
    assert float(summary["matched_flux_xy_max_incident_flux_p"]) > 0.0
    assert float(summary["counterpropagating_weighted_incident_total"]) == 0.0
    assert float(summary["counterpropagating_weighted_outgoing_total"]) == 0.0
    assert float(summary["counterpropagating_direct_incident_total"]) > 0.0
    assert float(summary["counterpropagating_direct_outgoing_total"]) > 0.0
    assert float(summary["counterpropagating_direct_incident_p"]) > 0.0
    assert float(summary["counterpropagating_direct_outgoing_t"]) > 0.0
    assert float(summary["counterpropagating_matched_source_abs_sum"]) == 0.0
    assert float(summary["counterpropagating_matched_exchange_power"]) == 0.0
    assert float(summary["counterpropagating_direct_matched_source_abs_sum"]) > 0.0
    assert float(summary["counterpropagating_direct_matched_exchange_power"]) < 0.0
    assert (
        float(summary["counterpropagating_propagated_launch_weighted_incident_total"])
        == 0.0
    )
    assert (
        float(summary["counterpropagating_propagated_launch_weighted_outgoing_total"])
        == 0.0
    )
    assert (
        float(summary["counterpropagating_propagated_launch_direct_incident_total"])
        > 0.0
    )
    assert (
        float(summary["counterpropagating_propagated_launch_direct_outgoing_total"])
        > 0.0
    )
    assert (
        float(summary["counterpropagating_propagated_launch_direct_incident_p"])
        > 0.0
    )
    assert (
        float(summary["counterpropagating_propagated_launch_direct_outgoing_t"])
        > 0.0
    )
    assert float(summary["counterpropagating_propagated_direct_final_energy"]) < float(
        summary["counterpropagating_propagated_matched_final_energy"]
    )
    assert float(
        summary["counterpropagating_propagated_direct_boundary_momentum"]
    ) < float(summary["counterpropagating_propagated_matched_boundary_momentum"])


def test_aethergrid_supports_classic_source_and_detector_hooks():
    grid = fdtd.AetherGrid(shape=(5, 5, 5))
    grid[2, 2, 2] = fdtd.PointSource(amplitude=1.0, phase_shift=pi / 2)
    grid[1:3, 1:3, 1:3] = fdtd.BlockDetector(name="detector")

    grid.step()

    assert float(grid.E[2, 2, 2, 2]) > 0.0
    assert len(grid.detector.E) == 1
    assert len(grid.detector.H) == 1


def test_aethergrid_supports_native_angular_detector_hook():
    grid = fdtd.AetherGrid(
        shape=(4, 4, 4),
        grid_spacing=1.0e6,
        native_angular_transport=True,
        native_angular_update_clocks=True,
    )
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid[0, :, :] = fdtd.AetherAngularSpongeBoundary(
        damping_t=0.25,
        damping_p=0.15,
    )
    grid[0:0, 0:0, 0:0] = fdtd.AetherNativeAngularDetector(
        name="native_detector",
        record_energy=True,
    )

    grid.step()

    assert len(grid.native_detector.readings["angular_momentum_t"]) == 1
    detected_momentum_t = np.asarray(
        grid.native_detector.readings["angular_momentum_t"][0]
    )
    detected_momentum_p = np.asarray(
        grid.native_detector.readings["angular_momentum_p"][0]
    )
    detected_omega_t = np.asarray(grid.native_detector.readings["omega_t"][0])
    detected_omega_p = np.asarray(grid.native_detector.readings["omega_p"][0])
    expected_momentum_t = float(grid.angular_momentum_t[0, 0, 0, 0])
    expected_momentum_p = float(grid.angular_momentum_p[0, 0, 0, 0])

    np.testing.assert_allclose(detected_momentum_t[0, 0, 0, 0], expected_momentum_t)
    np.testing.assert_allclose(detected_momentum_p[0, 0, 0, 0], expected_momentum_p)
    np.testing.assert_allclose(detected_omega_t[0, 0, 0, 0], expected_momentum_t / 2.0)
    np.testing.assert_allclose(detected_omega_p[0, 0, 0, 0], expected_momentum_p / 4.0)
    expected_energy = 0.5 * expected_momentum_t**2 / 2.0
    expected_energy += 0.5 * expected_momentum_p**2 / 4.0
    np.testing.assert_allclose(float(grid.native_detector.energy[0]), expected_energy)


def test_aethergrid_save_data_accepts_native_angular_detector_keys(tmp_path):
    grid = fdtd.AetherGrid(shape=(4, 4, 4), native_angular_transport=True)
    grid.folder = str(tmp_path)
    grid.angular_inertia_t[:, :, :, 0] = 2.0
    grid.angular_inertia_p[:, :, :, 0] = 4.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid[0:0, 0:0, 0:0] = fdtd.AetherNativeAngularDetector(
        name="native_detector",
        fields=("angular_momentum_t", "angular_momentum_p"),
        record_energy=True,
    )

    grid.step()
    grid.save_data()

    saved = np.load(tmp_path / "detector_readings.npz")
    assert "native_detector (angular_momentum_t)" in saved.files
    assert "native_detector (angular_momentum_p)" in saved.files
    assert "native_detector (native_angular_energy)" in saved.files


def test_aether_native_angular_detector_energy_requires_positive_inertia():
    grid = fdtd.AetherGrid(shape=(2, 2, 2))
    grid.angular_inertia_t[:, :, :, 0] = 1.0
    grid.angular_inertia_p[:, :, :, 0] = 0.0
    grid.angular_momentum_t[:, :, :, 0] = 0.75
    grid.angular_momentum_p[:, :, :, 0] = 1.25
    grid[0:0, 0:0, 0:0] = fdtd.AetherNativeAngularDetector(
        name="native_detector",
        record_energy=True,
    )

    try:
        grid.step()
    except ValueError as exc:
        assert "inertia" in str(exc)
    else:
        raise AssertionError("expected positive-inertia validation")

    assert grid.time_steps_passed == 0
    assert len(grid.native_detector.energy) == 0


def test_aethergrid_reset_clears_derived_fields():
    grid = fdtd.AetherGrid(shape=(5, 5, 5))
    grid[2, 2, 2] = fdtd.AetherPointSource(amplitude=1.0, phase_shift=pi / 2)

    grid.step()
    grid.reset()

    assert not np.any(np.asarray(grid.v))
    assert not np.any(np.asarray(grid.yank))
    assert grid.time_steps_passed == 0
