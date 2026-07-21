"""Regression tests for the migrated STPT PotentialGrid reference path."""

import numpy as np

from fdtd import PotentialGrid
from fdtd.operators import curl_face_to_edge, divergence, gradient


def test_potentialgrid_tracks_local_density_and_face_collocation():
    grid = PotentialGrid(4, 4, 4, dx=1.0, viscosity=1.0, background_density=2.0)

    grid.density[1, 1, 1] = 4.0
    grid.update_medium_fields()

    assert grid.reference_circulation == 0.5
    assert grid.local_circulation[1, 1, 1] == 0.25
    assert grid.face_density_x.shape == grid.vx.shape
    assert grid.face_density_y.shape == grid.vy.shape
    assert grid.face_density_z.shape == grid.vz.shape
    assert grid.face_circulation_x.shape == grid.vx.shape


def test_potentialgrid_step_runs_with_local_density_wave_solvers():
    grid = PotentialGrid(4, 4, 4, dx=1.0, viscosity=1.0, background_density=2.0)
    grid.density[1, 1, 1] = 4.0
    grid.vx[2, 2, 2] = 1.0

    grid.step()

    assert grid.step_count == 1
    assert grid.second_sound.face_density_x.shape == grid.second_sound.ax.shape
    assert np.isfinite(np.asarray(grid.density)).all()
    assert np.isfinite(np.asarray(grid.ax)).all()


def test_legacy_staggered_operator_dispatch_remains_available():
    scalar = np.zeros((4, 4, 4), dtype=float)
    scalar[1, 1, 1] = 1.0

    gx, gy, gz = gradient(scalar, dx=1.0)
    div = divergence(gx, gy, gz, dx=1.0)
    wx, wy, wz = curl_face_to_edge(gx, gy, gz, 1.0)

    assert gx.shape == (5, 4, 4)
    assert gy.shape == (4, 5, 4)
    assert gz.shape == (4, 4, 5)
    assert div.shape == scalar.shape
    assert wx.shape == (4, 5, 5)
    assert wy.shape == (5, 4, 5)
    assert wz.shape == (5, 5, 4)
