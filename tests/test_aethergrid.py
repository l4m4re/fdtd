"""Smoke tests for the experimental AetherGrid."""

import fdtd


def test_aethergrid_is_exported():
    assert hasattr(fdtd, "AetherGrid")


def test_aethergrid_construct_and_step():
    grid = fdtd.AetherGrid(shape=(4, 4, 4))
    grid.step()
    assert grid.time_steps_passed == 1
    assert grid.v.shape == (4, 4, 4, 3)
    assert grid.a.shape == (4, 4, 4, 3)
    assert grid.j.shape == (4, 4, 4, 3)
    assert grid.yank.shape == (4, 4, 4, 3)
