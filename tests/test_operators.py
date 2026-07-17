"""Tests for explicit aether bridge operators."""

import numpy as np

from fdtd.operators import (
    angular_clock_eigenvalue,
    angular_to_linear_bridge,
    apply_angular_metric,
    curl_edge_to_face,
    curl_face_to_edge,
    linear_to_angular_bridge,
)


def test_linear_to_angular_bridge_matches_legacy_curl_without_metric():
    field = np.zeros((4, 4, 4, 3), dtype=float)
    field[1, 1, 1, 0] = 2.0
    field[1, 2, 1, 2] = -1.0

    expected = curl_edge_to_face(field)
    actual = linear_to_angular_bridge(field)

    np.testing.assert_allclose(np.asarray(actual), np.asarray(expected))


def test_angular_to_linear_bridge_matches_legacy_curl_without_metric():
    field = np.zeros((4, 4, 4, 3), dtype=float)
    field[1, 1, 1, 1] = 3.0
    field[2, 1, 1, 2] = -0.5

    expected = curl_face_to_edge(field)
    actual = angular_to_linear_bridge(field)

    np.testing.assert_allclose(np.asarray(actual), np.asarray(expected))


def test_apply_angular_metric_weights_field_componentwise():
    field = np.ones((2, 2, 2, 3), dtype=float)
    metric_length = np.full((2, 2, 2, 1), 2.5, dtype=float)

    weighted = apply_angular_metric(field, metric_length)

    np.testing.assert_allclose(np.asarray(weighted), 2.5 * np.ones((2, 2, 2, 3)))


def test_angular_clock_eigenvalue_matches_small_step_limit():
    omega_t = np.array([1.25, 2.0, 0.5])
    omega_p = np.array([0.75, 0.25, 1.5])
    gamma = np.array([0.5, 1.0, 0.125])
    delta = 1e-5

    actual = angular_clock_eigenvalue(omega_t, omega_p, gamma, delta)
    expected = -(omega_t**2) - (omega_p**2) + gamma**2

    np.testing.assert_allclose(np.asarray(actual), expected, rtol=1e-10, atol=1e-10)


def test_angular_clock_eigenvalue_uses_exact_finite_step_formula():
    omega_t = 1.2
    omega_p = 0.8
    gamma = 0.4
    delta = 0.25

    actual = angular_clock_eigenvalue(omega_t, omega_p, gamma, delta)
    expected = (
        -4.0 / delta**2 * np.sin(delta * omega_t / 2.0) ** 2
        -4.0 / delta**2 * np.sin(delta * omega_p / 2.0) ** 2
        +4.0 / delta**2 * np.sinh(delta * gamma / 2.0) ** 2
    )

    np.testing.assert_allclose(np.asarray(actual), expected)


def test_angular_clock_eigenvalue_treats_gamma_as_independent_channel():
    delta = 0.25

    actual = angular_clock_eigenvalue(omega_t=0.0, omega_p=0.0, gamma=0.4, delta=delta)
    expected = 4.0 / delta**2 * np.sinh(delta * 0.4 / 2.0) ** 2

    np.testing.assert_allclose(np.asarray(actual), expected)


def test_angular_clock_eigenvalue_converges_to_small_step_limit():
    omega_t = 1.2
    omega_p = 0.9
    gamma = 0.35
    expected = -(omega_t**2) - (omega_p**2) + gamma**2

    coarse = float(angular_clock_eigenvalue(omega_t, omega_p, gamma, delta=0.1))
    fine = float(angular_clock_eigenvalue(omega_t, omega_p, gamma, delta=0.05))

    coarse_error = abs(coarse - expected)
    fine_error = abs(fine - expected)

    assert fine_error < coarse_error
    assert fine_error / coarse_error < 0.3


def test_angular_clock_eigenvalue_rejects_zero_step():
    try:
        angular_clock_eigenvalue(1.0, 1.0, 1.0, 0.0)
    except ValueError as exc:
        assert "delta" in str(exc)
    else:
        raise AssertionError("expected ValueError for zero delta")
