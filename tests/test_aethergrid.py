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


def test_aethergrid_step_does_not_promote_cartesian_omega_to_native_clocks():
    grid = fdtd.AetherGrid(shape=(5, 5, 5))
    grid[2, 2, 2] = fdtd.AetherPointSource(amplitude=1.0, phase_shift=pi / 2)

    grid.step()

    assert np.any(np.asarray(grid.angular_omega) != 0.0)
    assert not np.any(np.asarray(grid.omega_t))
    assert not np.any(np.asarray(grid.omega_p))
    assert not np.any(np.asarray(grid.angular_gamma))


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
