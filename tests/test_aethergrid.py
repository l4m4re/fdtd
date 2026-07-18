"""Tests for the experimental aether integration path."""

from math import pi

import numpy as np

import fdtd
from fdtd.operators import div


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
    assert grid.native_angular_momentum_rhs_t.shape == (5, 5, 5, 1)
    assert grid.native_angular_momentum_rhs_p.shape == (5, 5, 5, 1)
    assert grid.native_angular_momentum_candidate_t.shape == (5, 5, 5, 1)
    assert grid.native_angular_momentum_candidate_p.shape == (5, 5, 5, 1)


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
    assert not np.any(np.asarray(grid.native_angular_momentum_rhs_t))
    assert not np.any(np.asarray(grid.native_angular_momentum_rhs_p))
    assert not np.any(np.asarray(grid.native_angular_momentum_candidate_t))
    assert not np.any(np.asarray(grid.native_angular_momentum_candidate_p))


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
    grid.angular_inertia_t[1, 1, 1, 0] = 4.0
    grid.angular_inertia_p[1, 1, 1, 0] = 5.0
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
    assert not np.any(np.asarray(grid.native_angular_momentum_rhs_t))
    assert not np.any(np.asarray(grid.native_angular_momentum_rhs_p))
    assert not np.any(np.asarray(grid.native_angular_momentum_candidate_t))
    assert not np.any(np.asarray(grid.native_angular_momentum_candidate_p))


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


def test_aethergrid_supports_classic_source_and_detector_hooks():
    grid = fdtd.AetherGrid(shape=(5, 5, 5))
    grid[2, 2, 2] = fdtd.PointSource(amplitude=1.0, phase_shift=pi / 2)
    grid[1:3, 1:3, 1:3] = fdtd.BlockDetector(name="detector")

    grid.step()

    assert float(grid.E[2, 2, 2, 2]) > 0.0
    assert len(grid.detector.E) == 1
    assert len(grid.detector.H) == 1


def test_aethergrid_reset_clears_derived_fields():
    grid = fdtd.AetherGrid(shape=(5, 5, 5))
    grid[2, 2, 2] = fdtd.AetherPointSource(amplitude=1.0, phase_shift=pi / 2)

    grid.step()
    grid.reset()

    assert not np.any(np.asarray(grid.v))
    assert not np.any(np.asarray(grid.yank))
    assert grid.time_steps_passed == 0
