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

Native angular staging
----------------------

The current bridge exposes a second, passive angular staging layer alongside
the Cartesian ``angular_omega`` field. This layer is meant to make the intended
native angular geometry testable without silently changing the timestep update.

The staged native fields are:

- ``omega_t`` and ``omega_p`` for two angular clock channels;
- ``angular_gamma`` for a hyperbolic or envelope benchmark channel;
- ``theta_t``, ``theta_p``, and ``angular_chi`` for the corresponding
  dimensionless timestep increments;
- ``angular_clock_lambda`` for the local angular-clock eigenvalue benchmark;
- ``angular_e_t`` and ``angular_e_p`` for the current local angular frame;
- ``ell_t`` and ``ell_p`` for lever-arm or metric lengths;
- ``angular_inertia_t`` and ``angular_inertia_p`` for provisional native
  moment-density channels;
- ``angular_momentum_t`` and ``angular_momentum_p`` for native angular
  momentum channels; and
- ``angular_torque_t`` and ``angular_torque_p`` for staged native torque
  channels.

The residual path also stores ``native_angular_source_t`` and
``native_angular_source_p`` as explicit source or boundary-exchange terms.
These arrays are zero by default and are not inferred from ordinary Maxwell
boundary behavior.

For passive transport experiments, ``native_angular_momentum_rhs_t`` and
``native_angular_momentum_rhs_p`` store the candidate right-hand side
``S_L - ell*div(tau_native)``. ``native_angular_momentum_candidate_t`` and
``native_angular_momentum_candidate_p`` store one predicted next momentum state
without promoting it into ``angular_momentum_t`` or ``angular_momentum_p``.

The helper methods that populate these fields are deliberately explicit:

- ``evaluate_angular_clock_benchmark()`` evaluates
  ``omega_t``, ``omega_p``, and ``angular_gamma`` without changing the bridge
  dynamics;
- ``update_native_angular_inertia()`` applies the provisional closure
  ``I = rho * ell**2``;
- ``update_native_angular_momentum()`` applies ``L = I * omega``;
- ``update_native_angular_torque()`` evaluates a finite-difference
  ``tau = dL/dt`` from a previous momentum state; and
- ``project_native_angular_torque()`` maps the two native torque channels into
  ``native_angular_tau`` through the local frame vectors;
- ``evaluate_native_angular_momentum_rhs()`` evaluates the passive transport
  right-hand side for ``dL/dt``; and
- ``predict_native_angular_momentum_step()`` computes a candidate next native
  angular momentum state without applying it.

None of these helpers are called by ``step()``. They are diagnostics and
staging points for the next architecture pass, not a completed angular update
law.

Residual diagnostics
--------------------

The native torque projection can also be inspected spatially:

- ``evaluate_native_angular_torque_divergence()`` computes
  ``div(native_angular_tau)``;
- ``evaluate_native_angular_metric_divergence()`` multiplies that divergence
  by ``ell_t`` and ``ell_p`` so it can be compared with the staged torque
  channels; and
- ``evaluate_native_angular_transport_residual()`` evaluates the passive
  candidate

``R_L = dL/dt + ell*div(tau_native) - S_L``.

The residual helper treats ``angular_torque_t`` and ``angular_torque_p`` as the
``dL/dt`` terms. Any external or boundary exchange must be supplied explicitly
through the optional source arrays. The current implementation does not infer
boundary exchange, advance angular momentum, or feed this residual back into
``linear_a``.

Equivalently, the passive transport RHS is
``dL/dt = S_L - ell*div(tau_native)``, and the residual is the difference
between the staged torque and that RHS. This is useful for bounded diagnostic
tests because a caller can compare a staged torque law with the candidate
transport law before any update rule is promoted.

The first bounded predictor benchmark is the balanced-source no-drift case in
``test_aethergrid_passive_momentum_predictor_balanced_64_step_benchmark``. It
uses 64 candidate-only iterations with explicit sources matching the
metric-weighted torque divergence. The candidate momentum fields must stay
finite and equal to the initial native angular momentum. This benchmark does
not promote the predictor into ``step()`` and does not define physical boundary
exchange.

The first unbalanced accounting benchmark is
``test_aethergrid_passive_momentum_predictor_damping_64_step_benchmark``. It
adds explicit local damping sources and requires the passive candidate momentum
to follow the exact discrete decay. This is still only source accounting; it is
not a sponge layer, absorbing boundary, or physical boundary law.

The first sponge-style accounting benchmark is
``test_aethergrid_passive_momentum_predictor_sponge_64_step_benchmark``. It
uses explicit boundary-layer damping masks in the source arrays so only masked
candidate momentum cells decay. It still does not call boundary hooks or define
absorbing or reflective boundary physics.

For opt-in experiments, scene elements may expose
``native_angular_source_terms()`` and return ``(source_t, source_p)`` arrays.
``collect_native_angular_source_terms()`` sums those explicit terms into
``native_angular_source_t`` and ``native_angular_source_p``. This collector is
not called by ``step()`` or by the residual helper; callers must pass the
collected arrays to ``evaluate_native_angular_transport_residual()`` when that
is the intended test.

Caveat
------

The longer-term STPT direction in this repository is stronger than a simple
reuse of Maxwell staggering. The working hypothesis is that linear and angular
branches may need distinct but coupled discrete sectors. ``AetherGrid`` should
therefore be understood as an integration bridge: useful for scene
construction, comparison against Maxwell examples, and backend reuse, but not
yet the finished geometric simulator.

The next unresolved design step is still physical boundary semantics for the
native angular sector. The collector only defines how explicit exchange terms
can enter the diagnostic residual. A future update rule must state where
angular boundary data live, which of clocks, momentum, torque, projected
torque, or residuals a boundary may modify, how exchange enters the explicit
source terms, and what bounded-growth test accepts the result.

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
