"""Tests for explicit aether bridge operators."""

import numpy as np

from fdtd.operators import (
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
