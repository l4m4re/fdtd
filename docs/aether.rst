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
``native_angular_linear_response`` stores a passive projection of staged native
angular stress back toward the linear-sector layout without writing
``angular_A`` or ``linear_a``.
The charge-observable protocol is staged as passive bookkeeping too:
``linear_charge_flux_candidate`` stores the normalized local linear integrand,
``native_angular_charge_flux_t`` and ``native_angular_charge_flux_p`` store the
two native angular loop-channel candidates, and
``native_angular_charge_reduction`` stores one selected angular reduction for
comparison.

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
  angular momentum state without applying it; and
- ``advance_native_angular_momentum_transport()`` promotes that candidate into
  ``angular_momentum_t`` and ``angular_momentum_p`` only when called
  explicitly;
- ``evaluate_native_angular_linear_response()`` projects ``native_angular_tau``
  through the current angular-to-linear bridge as a coupling diagnostic without
  applying feedback;
- ``evaluate_linear_charge_flux_candidate()`` evaluates
  ``s_l*delta_t*eta*v*dA/A0`` as a local linear candidate;
- ``evaluate_native_angular_charge_flux_candidates()`` evaluates the signed
  ``delta_t*eta*ell*omega`` two-clock angular candidates; and
- ``reduce_native_angular_charge_flux()`` compares the additive and
  geometric-mean angular reductions without choosing which is physical.

None of these helpers are called by ``step()``. They are diagnostics and
staging points for the next architecture pass, not a completed production
angular update law.

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
through the optional source arrays. The residual and predictor helpers do not
infer boundary exchange or feed this residual back into ``linear_a``.
Previously collected source buffers are not read implicitly by the residual or
predictor helpers.
The current spatial transport route is guarded by analytic stencil checks:
uniform native torque must have zero divergence, and a linearly varying native
torque field must reproduce the exact finite-difference pattern used by the
collocated prototype before the residual is balanced.
The first bounded propagated-transport check is still candidate-only: it
derives staged native torque from successive momentum candidates, projects that
torque, evaluates the metric divergence, and feeds only the next candidate
state for a fixed number of iterations. It must remain finite and uncoupled
from ``step()`` and ``linear_a``.
``advance_native_angular_momentum_transport()`` is the corresponding opt-in
transport advance. It promotes the predicted native angular momentum candidate
through ``L_next = L + dt*(S_L - ell*div(tau_native))``, can collect explicit
opt-in source and boundary hooks when requested, and can update the native
clocks from positive inertia when requested. It is still not called by
``step()`` and does not update the legacy bridge fields.
The first linear/angular coupling gate is similarly passive:
``evaluate_native_angular_linear_response()`` must match the current
``angular_to_linear_bridge`` projection while leaving ``angular_A`` and
``linear_a`` unchanged.

Equivalently, the passive transport RHS is
``dL/dt = S_L - ell*div(tau_native)``, and the residual is the difference
between the staged torque and that RHS. This is useful for bounded diagnostic
tests because a caller can compare a staged torque law with the candidate
transport law before the update is allowed into the production timestep.

The first bounded predictor benchmark is the balanced-source no-drift case in
``test_aethergrid_passive_momentum_predictor_balanced_64_step_benchmark``. It
uses 64 candidate-only iterations with explicit sources matching the
metric-weighted torque divergence. The candidate momentum fields must stay
finite and equal to the initial native angular momentum. This benchmark does
not promote the predictor into ``step()`` and does not define physical boundary
exchange.
``test_aethergrid_native_angular_transport_advance_balanced_64_step_benchmark``
applies the same no-drift gate to
``advance_native_angular_momentum_transport()``. The helper may promote native
angular momentum, but the promoted state must remain finite, equal to the
initial state under balanced sources, and uncoupled from ``step()``,
``angular_tau``, ``angular_A``, and ``linear_a``.

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

The first opt-in boundary hook is ``AetherAngularSpongeBoundary``. It exposes
``native_angular_source_terms()`` and contributes ``-damping*L`` exchange terms
on its registered grid slice when callers explicitly collect native angular
source terms. It is a passive damping hook, not a Maxwell PML, absorber, or
reflector.
``AetherAngularNoExchangeBoundary`` is the matching passive closed-boundary
baseline: it exposes the same hook but returns zero exchange arrays. It is not
a reflection law; it only makes "boundary present, no native-angular exchange"
explicit in the source-term contract.
``AetherAngularReflectiveBoundary`` adds a passive mirror-relaxation hook: it
uses the adjacent interior cell as the target state and contributes
``rate*(sign*L_mirror - L_boundary)`` on the registered boundary slice. This is
valid only when registered on one outer grid face; non-outer slices are
rejected. This is still diagnostic source accounting, not a production
wave-reflection boundary.
``AetherGrid.evaluate_native_angular_kinetic_energy()`` adds the matching
quadratic bookkeeping diagnostic ``0.5*L**2/I`` for the native angular
momentum channels. It is used to check candidate-only source accounting, not to
define a Hamiltonian, detector observable, or production energy update.
The sponge-boundary energy benchmark uses the same diagnostic to require
monotone decay under local ``-damping*L`` exchange on boundary slices. That is
still damping bookkeeping; it is not yet an absorbing boundary proof.

For opt-in experiments, scene elements may expose
``native_angular_source_terms()`` and return ``(source_t, source_p)`` arrays.
``collect_native_angular_source_terms()`` sums those explicit terms into
``native_angular_source_t`` and ``native_angular_source_p``. This collector is
not called by ``step()`` or by the residual helper; callers must pass the
collected arrays to ``evaluate_native_angular_transport_residual()`` when that
is the intended test. The same explicit-passing rule applies to
``predict_native_angular_momentum_step()``. The collector resets both source
buffers on every call, so ``include_sources`` and ``include_boundaries`` select
a fresh accounting view instead of accumulating previous collections.
The current native-angular boundary contract is intentionally narrow: boundary
hooks may read native angular momentum and return source arrays, but they may
not mutate clocks, momentum, torque fields, residuals, candidates, or
linear-sector state.

Charge-flux candidates
----------------------

The current bridge exposes the charge-flux protocol as explicit candidate
bookkeeping, not as an experimental fit to measured charge.
``evaluate_linear_charge_flux_candidate()`` computes the local linear candidate
``s_l*delta_t*eta*v*dA/A0`` and stores it in
``linear_charge_flux_candidate``. The caller must supply the control-surface
measure and reference area when using anything other than the default local
unit measure.

``evaluate_native_angular_charge_flux_candidates()`` computes the native
angular channel candidates
``s_t*delta_t*eta*ell_t*omega_t*dtheta_t/(2*pi)/N_a`` and
``s_p*delta_t*eta*ell_p*omega_p*dtheta_p/(2*pi)/N_a`` and stores them in
``native_angular_charge_flux_t`` and ``native_angular_charge_flux_p``.
``reduce_native_angular_charge_flux()`` then stores either ``q_t+q_p`` or
``sqrt(q_t*q_p)`` in ``native_angular_charge_reduction``. The geometric-mean
candidate requires the chosen signs to make the channel product non-negative.

These helpers make ``A0``, ``delta_t``, signs, loop measures, angular
normalization, and reduction mode explicit at the call site. They do not
derive the elementary charge, select a physical charge polarity convention, or
feed any charge observable back into the timestep update.

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
