"""Tests for the experimental aether integration path."""

from math import pi

import numpy as np

import fdtd


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
