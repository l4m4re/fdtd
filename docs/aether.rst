Experimental Aether Path
========================

The public :mod:`fdtd` package is a Maxwell FDTD simulator first. The
experimental :class:`fdtd.AetherGrid` path keeps the same package rhythm while
exploring a different theory:

- the scene is still built with ``grid[...] = ...``;
- sources, boundaries, detectors, and objects still register on the grid;
- the same backend abstraction is reused for NumPy, Torch, and CUDA backends;
- and the existing Maxwell examples remain the natural comparison baseline.

What changes is the state being evolved. ``AetherGrid`` currently advances a
velocity-like field ``v`` and derives scalar, rotational, electromagnetic,
acceleration, and jerk-like fields from it. The present implementation is still
an educational prototype and should be read as an explicit stepping stone
toward a more geometric STPT simulator rather than as a final theory claim.

Current scope
-------------

The aether path is intentionally opt-in:

- classic :class:`fdtd.Grid` behavior is left untouched;
- aether-native sources can inject directly into ``v`` through
  :class:`fdtd.AetherPointSource` and :class:`fdtd.AetherLineSource`;
- Maxwell-style sources and detectors can still be reused for side-by-side
  comparison because ``AetherGrid`` exposes ``E`` and ``H``-like derived
  fields.

Operator picture
----------------

The current implementation reuses the package's discrete operator patterns in a
more neutral language:

- ``curl_edge_to_face`` maps circulation from edge-like slots to face-like
  slots;
- ``curl_face_to_edge`` maps the reverse direction;
- ``linear_to_angular_bridge`` names the current temporary map from the linear
  translational sector into the angular sector;
- ``angular_to_linear_bridge`` names the reverse temporary coupling;
- ``apply_angular_metric`` marks where the angular sector is expected to carry
  an explicit meter-valued length or lever-arm factor;
- ``div`` collapses an edge-like vector field to a scalar slot; and
- ``grad`` lifts a scalar slot back to a vector field.

This keeps the code close to the original educational package style while the
theory is still being sorted out.

``AetherGrid`` now also keeps an explicit ``angular_metric_length`` field in
its internal state. The current runnable bridge does not yet feed that metric
back into the dynamics, but the placeholder is there so the dimensional
separation between linear transport and angular geometry remains visible in the
code.

Caveat
------

The longer-term STPT direction in this repository is stronger than a simple
reuse of Maxwell staggering. The working hypothesis is that linear and angular
branches may need distinct but coupled discrete sectors. ``AetherGrid`` should
therefore be understood as an integration bridge: useful for scene
construction, comparison against Maxwell examples, and backend reuse, but not
yet the finished geometric simulator.

Quick-start adaptation
----------------------

There is now a small runnable adaptation of the original quick-start scene in
``examples/aether_quick_start.py``. It keeps the same source placement,
detector placement, and overall scene-construction style, but switches the
native source to :class:`fdtd.AetherLineSource`.

Run it from the package root with for example::

    python examples/aether_quick_start.py --maxwell-steps 40 --aether-steps 8

The default aether run uses fewer steps and a much smaller source amplitude
than the Maxwell run. That is intentional: the current prototype is useful for
studying how the scene is wired and where the dynamics become unstable, but it
is not yet numerically mature enough to treat the stock Maxwell quick-start
parameters as physically neutral.
