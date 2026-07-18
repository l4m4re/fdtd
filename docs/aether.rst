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
- :class:`fdtd.AetherNativeAngularPointSource` can provide explicit
  native-angular ``S_t/S_p`` source terms for opt-in transport experiments
  without injecting into Maxwell ``E/H`` fields;
- Maxwell-style sources and detectors can still be reused for side-by-side
  comparison because ``AetherGrid`` exposes ``E`` and ``H``-like derived
  fields; and
- :class:`fdtd.AetherNativeAngularDetector` can record opt-in native-angular
  momentum, clock, source, power, and local energy diagnostics without using
  the Maxwell ``E/H`` detector contract.

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
``native_angular_exchange_power_t`` and
``native_angular_exchange_power_p`` store the companion source-power
diagnostic ``S_L*L/I`` for classifying explicit exchange as injecting,
neutral, or dissipative.
``native_angular_transport_power_t`` and
``native_angular_transport_power_p`` store the full transport-RHS diagnostic
``(S_L - ell*div(tau_native))*L/I``.
``native_angular_boundary_direct_*_flux_t`` and
``native_angular_boundary_direct_*_flux_p`` store the direct native-channel
boundary comparator, where t and p are projected and sign-split before their
normal fluxes are summed.

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
- ``configure_native_angular_frame()`` sets a uniform or per-cell
  ``angular_e_t/p`` frame and validates nonzero, orthogonal axes;
- ``project_native_angular_torque()`` maps the two native torque channels into
  ``native_angular_tau`` through the local frame vectors;
- ``evaluate_native_angular_exchange_power()`` evaluates ``S_L*L/I`` for
  explicit source or boundary exchange;
- ``evaluate_native_angular_boundary_flux()`` projects ``native_angular_tau``
  onto a selected outer-face normal and splits it into incident and outgoing
  diagnostic flux;
- ``collect_native_angular_boundary_fluxes()`` applies that flux diagnostic to
  registered boundary contracts that advertise outer-face metadata;
- ``evaluate_native_angular_boundary_channel_flux()`` and
  ``collect_native_angular_boundary_channel_fluxes()`` weight the same
  incident/outgoing flux by the local ``angular_e_t/p`` normal projection;
- ``evaluate_native_angular_boundary_direct_channel_flux()`` and
  ``collect_native_angular_boundary_direct_channel_fluxes()`` project and
  sign-split ``tau_t`` and ``tau_p`` separately before summing them, giving a
  comparator for the matched-flux candidate;
- ``evaluate_native_angular_momentum_rhs()`` evaluates the passive transport
  right-hand side for ``dL/dt``; and
- ``evaluate_native_angular_transport_power()`` evaluates the corresponding
  full-RHS power diagnostic;
- ``predict_native_angular_momentum_step()`` computes a candidate next native
  angular momentum state without applying it; and
- ``advance_native_angular_momentum_transport()`` promotes that candidate into
  ``angular_momentum_t`` and ``angular_momentum_p`` only when called
  explicitly;
- ``evaluate_native_angular_linear_response()`` projects ``native_angular_tau``
  through the current angular-to-linear bridge as a coupling diagnostic without
  applying feedback;
- ``apply_native_angular_linear_response()`` applies that response to
  ``linear_a`` only when called explicitly, with an ``add`` or ``replace``
  mode chosen by the caller;
- ``evaluate_linear_charge_flux_candidate()`` evaluates
  ``s_l*delta_t*eta*v*dA/A0`` as a local linear candidate;
- ``evaluate_native_angular_charge_flux_candidates()`` evaluates the signed
  ``delta_t*eta*ell*omega`` two-clock angular candidates; and
- ``reduce_native_angular_charge_flux()`` compares the additive and
  geometric-mean angular reductions without choosing which is physical.

By default, these helpers are not called by ``step()``. When ``AetherGrid`` is
constructed with ``native_angular_transport=True`` or a non-``None``
``native_angular_feedback_mode``, ``update()`` runs the corresponding
native-angular transport or feedback hook inside the normal timestep lifecycle.
That path is opt-in and still not a completed physical angular update law.

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
``test_aethergrid_step_dynamic_native_torque_transport_loop_is_bounded`` is
the first flagged ``step()`` version of that harness. It derives native torque
from the previous promoted momentum state before each step, lets the opt-in
transport lifecycle promote the next state, verifies clock resync from positive
inertia, records native detector samples, and keeps the linear velocity and
explicit source buffers quiet.
``advance_native_angular_momentum_transport()`` is the corresponding opt-in
transport advance. It promotes the predicted native angular momentum candidate
through ``L_next = L + dt*(S_L - ell*div(tau_native))``, can collect explicit
opt-in source and boundary hooks when requested, and can update the native
clocks from positive inertia when requested. It is called by ``step()`` only
when ``native_angular_transport=True`` and does not update the legacy bridge
fields.
The first linear/angular coupling gate starts as a passive diagnostic:
``evaluate_native_angular_linear_response()`` must match the current
``angular_to_linear_bridge`` projection while leaving ``angular_A`` and
``linear_a`` unchanged.
``apply_native_angular_linear_response()`` is the corresponding opt-in
feedback helper: it can add or replace ``linear_a`` with an explicit or freshly
projected response. It is called by ``step()`` only when
``native_angular_feedback_mode`` is ``"add"`` or ``"replace"``, and it still
does not write ``angular_A`` or define the physical feedback law.

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
``test_aethergrid_step_runs_opt_in_native_angular_transport`` then verifies the
same transport hook through ``step()`` with the constructor flag enabled, while
the default-step test verifies that native angular source collection and
momentum promotion remain disabled by default.
``test_aethergrid_step_native_angular_transport_no_drift_64_step`` is the
matching lifecycle bounded gate: 64 flagged ``step()`` calls with balanced
native source terms must leave native angular momentum finite and unchanged.
``test_aethergrid_step_native_angular_feedback_replace_64_step_is_bounded``
adds the feedback lifecycle gate: a fixed native response in ``replace`` mode
must remain finite, keep ``linear_a`` equal to the selected response, and drive
linear velocity predictably across 64 flagged ``step()`` calls.
``test_aethergrid_step_dynamic_native_torque_feedback_drives_bounded_linear_response``
adds the paired dynamic coupling gate: the same native torque transport loop is
run with and without replacement feedback, native momentum must match in both
runs, and only the feedback run may produce a finite nonzero linear velocity
response. This is still an opt-in bridge coupling check, not a physical
feedback law.

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
``test_aethergrid_transport_advance_sponge_boundary_reduces_quadratic_energy``
applies the sponge hook to the promoted opt-in transport helper and requires
the same local discrete damping curve plus monotone decrease of the quadratic
diagnostic.
``test_aethergrid_step_sponge_boundary_reduces_quadratic_energy_64_step``
checks the same damping curve through ordinary ``step()`` calls when
``native_angular_transport=True``.
``AetherAngularNoExchangeBoundary`` is the matching passive closed-boundary
baseline: it exposes the same hook but returns zero exchange arrays. It is not
a reflection law; it only makes "boundary present, no native-angular exchange"
explicit in the source-term contract.
``test_aethergrid_transport_advance_no_exchange_boundary_hook_64_step`` applies
that closed-boundary baseline to the promoted opt-in helper and requires a
balanced state to remain finite and unchanged.
``AetherAngularReflectiveBoundary`` adds a passive mirror-relaxation hook: it
uses the adjacent interior cell as the target state and contributes
``rate*(sign*L_mirror - L_boundary)`` on the registered boundary slice. This is
valid only when registered on one outer grid face; non-outer slices are
rejected. This is still diagnostic source accounting, not a production
wave-reflection boundary.
``AetherAngularMatchedFluxBoundary`` is the first candidate physical boundary
law. It reads the incident part of the projected ``native_angular_tau`` normal
flux, projects that boundary slot onto the local ``angular_e_t/p`` frame, and
returns source terms that oppose local native angular momentum. Its contract
sets ``physical_boundary_law=True`` and supplies the required boundary-slot,
channel-weighted flux-split, metric-frame, zero-incident-response,
energy-balance, falsification-comparator, and acceptance-test metadata.
``test_aether_angular_matched_flux_boundary_candidate_is_passive`` checks the
limited passivity claim: for positive inertia, the candidate's
``S_L*L/I`` exchange power is non-positive. This is still a candidate on the
collocated bridge, not a derived native t/p absorber. The current x-face
example deliberately exposes a falsification target: with the default
``angular_e_t=x`` and ``angular_e_p=y`` frame, x-face matched flux only has a
t-channel normal projection; p-channel absorption must be tested on y faces or
with a different angular frame. The swapped-frame variant uses
``configure_native_angular_frame()`` to exchange the x/y channel visibility on
the same x faces. In that swapped scene, final quadratic energy and
nonpositive exchange power remain the acceptance checks; boundary-band
momentum is reported as a frame-sensitive diagnostic rather than as an
invariant.
``test_aether_angular_matched_flux_boundary_candidate_no_incident_no_exchange``
adds the corresponding negative guard: if the registered matched-flux faces see
only outgoing projected normal flux, the candidate must return zero source
exchange and zero exchange power.
``test_aethergrid_direct_native_channel_flux_exposes_projected_cancellation``
and the boundary-flux example's ``counterpropagating_*`` output add the first
native t/p comparator case: two oblique orthogonal channels can cancel in the
summed projected normal flux while still carrying equal and opposite direct
channel flux. This is a falsification target for the current matched-flux
rule, not a replacement boundary law.
``AetherAngularDirectMatchedFluxBoundary`` is the matching direct-channel
candidate. It uses that comparator quantity as its source input: t and p are
projected onto the boundary normal, sign-split separately, and only incident
native-channel flux is converted into dissipative source exchange. Its current
role is to compare against ``AetherAngularMatchedFluxBoundary`` in blind-spot
cases, not to claim a finished absorber.
``test_aethergrid_transport_advance_reflective_boundary_conserves_quadratic_energy``
checks the same sign-flip energy bookkeeping through the promoted opt-in
transport helper.
``test_aethergrid_step_reflective_boundary_conserves_quadratic_energy_64_step``
then runs the same bookkeeping gate through ordinary ``step()`` calls when
``native_angular_transport=True``.
``test_aethergrid_reflective_boundary_can_change_energy_when_unmatched`` is the
matching negative guard: an unmatched boundary and mirror state can change the
quadratic diagnostic. The current reflective hook is therefore only
mirror-relaxation source accounting, not a general energy-conserving physical
reflection law.
``AetherGrid.evaluate_native_angular_kinetic_energy()`` adds the matching
quadratic bookkeeping diagnostic ``0.5*L**2/I`` for the native angular
momentum channels. It is used to check candidate-only source accounting, not to
define a Hamiltonian, detector observable, or production energy update.
The sponge-boundary energy benchmark uses the same diagnostic to require
monotone decay under local ``-damping*L`` exchange on boundary slices. That is
still damping bookkeeping; it is not yet an absorbing boundary proof.
``test_aether_angular_no_exchange_boundary_has_zero_exchange_power`` and
``test_aethergrid_step_sponge_boundary_matches_exchange_power_balance`` add
the first source-power gates: no-exchange boundaries must report zero power,
and sponge exchange through ordinary opt-in ``step()`` calls must match
``Delta E = dt*P + 0.5*dt**2*S**2/I`` for the staged quadratic diagnostic.
This is passivity bookkeeping, not a propagated-wave absorber proof.
``test_aethergrid_step_dynamic_transport_sponge_boundary_is_dissipative``
then compares the same staged dynamic transport loop with and without x-face
sponge hooks. The sponge run must have nonpositive exchange power and no
greater total quadratic energy or boundary-band momentum than the no-sponge
run. This is still a lifecycle passivity check, not a finished absorber.
``test_aethergrid_step_static_torque_transport_matches_power_balance`` applies
the same finite-step balance to the full transport RHS with nonzero static
native torque divergence. This is bounded transport bookkeeping, not a
physical wave-propagation claim.
``test_aethergrid_evaluates_native_angular_boundary_flux_split`` and
``test_aethergrid_collects_native_angular_boundary_fluxes_from_contracts`` add
the first incident/outgoing boundary diagnostic gates: the outward-normal
component of ``native_angular_tau`` is sign-split on selected outer faces and
can be evaluated through registered boundary contracts. This is a projected
bridge diagnostic, not a two-channel native angular flux law or a physical
absorber/reflector.
``test_aethergrid_evaluates_native_angular_boundary_channel_flux_split`` and
``test_aethergrid_collects_native_angular_boundary_channel_fluxes_from_contracts``
add the channel view of that diagnostic by weighting incident/outgoing flux
with the local ``angular_e_t/p`` normal projection.
``test_aethergrid_direct_native_channel_flux_matches_axis_aligned_projection``
and ``test_aethergrid_collects_direct_native_channel_fluxes_from_contracts``
add the direct native-channel comparator. It matches the weighted diagnostic in
axis-aligned cases and gives registered boundary-contract summaries without
collecting source terms.
``test_aethergrid_configures_uniform_native_angular_frame`` and
``test_aethergrid_configures_field_native_angular_frame`` keep that frame
configuration explicit and reject malformed, zero, or non-orthogonal inputs.
``test_aethergrid_step_native_angular_transport_updates_clocks_64_step`` checks
the optional clock-resync path: when ``native_angular_update_clocks=True`` and
native inertia is positive, ``omega_t`` and ``omega_p`` must track ``L/I``
through ordinary opt-in ``step()`` calls. Missing positive inertia raises
before the timestep counter advances.
``test_aethergrid_supports_native_angular_detector_hook`` adds the first
detector-facing native-angular observable gate: the detector registers through
``grid[...] =`` and samples after the opt-in native-angular update. The save
path is also generic over detector keys, so native-angular readings can be
stored without pretending they are Maxwell ``E`` or ``H`` samples.

For opt-in experiments, scene elements may expose
``native_angular_source_terms()`` and return ``(source_t, source_p)`` arrays.
``collect_native_angular_source_terms()`` sums those explicit terms into
``native_angular_source_t`` and ``native_angular_source_p``. This collector is
called by ``step()`` only when ``native_angular_transport=True`` and
``native_angular_collect_sources=True``. The residual and predictor helpers
still require callers to pass collected arrays explicitly when that is the
intended test. The collector resets both source buffers on every call, so
``include_sources`` and ``include_boundaries`` select a fresh accounting view
instead of accumulating previous collections.
``AetherNativeAngularPointSource`` is the first public source implementing
this contract. It is a source-term provider only: its Maxwell ``update_E`` and
``update_H`` hooks are no-ops, so any native-angular momentum change must come
through the opt-in transport collector.
The current native-angular boundary contract is intentionally narrow: boundary
hooks may read native angular momentum and return source arrays, but they may
not mutate clocks, momentum, torque fields, residuals, candidates, or
linear-sector state.
Every current native-angular boundary hook also exposes
``native_angular_boundary_contract()``. These contracts identify the hook as
``scope="source_accounting"``, list the returned source buffers, declare no
direct mutation targets, and set ``physical_boundary_law`` to ``False``.
Current hooks also report their registered outer face axis and side.
``AetherGrid.collect_native_angular_boundary_contracts()`` collects those hook
contracts at grid level and annotates them with boundary type and registered
name without collecting source buffers.
``AetherGrid.validate_native_angular_boundary_contracts()`` verifies that the
current hooks remain only source-accounting contracts, while
``require_physical_laws=True`` is reserved as the future acceptance gate for
real absorbers or reflectors. A future physical absorber or reflector should
expose a different, stronger contract rather than being inferred from these
source-accounting hooks.

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
source terms, how incident and outgoing angular momentum or torque flux are
separated, which metric lengths and orientation frames define the boundary
slots, and what propagated-wave bounded-growth test accepts the result. Until
that stronger contract exists, current hooks should continue to report
``physical_boundary_law=False`` and should fail validation when
``require_physical_laws=True``.

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

Native-angular observable scene
-------------------------------

``examples/aether_native_angular_observables.py`` is the smaller scene for the
native-angular staging path itself. It registers an
``AetherNativeAngularPointSource`` and an ``AetherNativeAngularDetector`` on a
tiny opt-in ``AetherGrid`` with positive native angular inertia.

Run it from the package root with for example::

    python examples/aether_native_angular_observables.py --steps 8

The expected acceptance is deliberately modest: source terms are collected by
the native-angular transport path, the detector records finite native
momentum, clock, source, power, and quadratic energy samples, and the linear
velocity field remains quiet. This is a reproducible observable benchmark for
the package-style API, not a propagated angular-wave or physical boundary
validation.

Dynamic native-angular transport scene
--------------------------------------

``examples/aether_native_angular_dynamic_transport.py`` is the runnable scene
for the current dynamic Phase 2 path. It derives native torque from the
previous promoted momentum layer before each step, projects that torque into
``native_angular_tau``, and then lets the opt-in ``AetherGrid.step()``
lifecycle promote the next native momentum state.

Run the most complete small example from the package root with::

    python examples/aether_native_angular_dynamic_transport.py --steps 16 --sponge --feedback --charge-observables

With those flags, the scene also enables the current x-face sponge
source-accounting hooks and replacement feedback from the native-angular
linear response into ``linear_v``. ``--charge-observables`` evaluates the
explicit local linear and native-angular charge-flux candidates after the run
using the normalizations supplied on the command line. The example also
reports the grid-level native-angular boundary contract count, whether all
reported contracts are source-accounting contracts, and whether any current
boundary hook claims to be a physical boundary law. It also reports the maximum
x-face incident, outgoing, and net projected boundary flux seen through the
registered boundary contracts during the run. The expected acceptance is
finite native momentum, finite linear velocity, nonzero staged transport, one
native detector sample per step, finite charge
candidates, visible boundary-flux diagnostics, and
``physical_boundary_laws: False`` for the current sponge setup. This is the
package-facing dynamic workflow benchmark for the bridge; it is not a finished
native angular lattice, absorber, reflection law, physical feedback
derivation, or independent charge derivation.

Boundary-flux acceptance scene
------------------------------

``examples/aether_native_angular_boundary_flux.py`` is the first dedicated
boundary-observability benchmark for Phase 2. It runs the same staged
native-angular transport setup with registered x-face
``AetherAngularNoExchangeBoundary`` hooks, x-face
``AetherAngularSpongeBoundary`` hooks, x-face
``AetherAngularMatchedFluxBoundary`` hooks, a swapped-frame x-face matched-flux
variant, and an xy-face no-exchange / matched-flux pair to expose both default
angular frame channels in one benchmark. All propagated runs use
contract-driven
``collect_native_angular_boundary_channel_fluxes()`` diagnostics.

Run it from the package root with for example::

    python examples/aether_native_angular_boundary_flux.py --steps 32

The expected acceptance is finite native momentum in all runs, two registered
boundary contracts per propagated run, ``physical_boundary_laws: False`` for
no-exchange and sponge, ``physical_boundary_laws: True`` for the matched-flux
candidate, visible incident/outgoing and matched-flux t/p channel diagnostics,
nonpositive sponge and matched-flux exchange power, and no greater final sponge
or matched-flux energy than the no-exchange baseline. Boundary-band momentum
is reported as a diagnostic, but it is not an invariant for the swapped-frame
scene. This is the current propagated-boundary acceptance instrumentation and
static candidate-comparison harness; neither candidate is yet a derived native
t/p angular absorber or reflector. With the current default frame, the x-face
matched-flux branch reports zero p-channel incident flux. The swapped-frame
x-face branch and the xy-face branch must both make p-channel incident flux
positive.
The same output also reports a static counterpropagating oblique-frame case:
``counterpropagating_weighted_incident_total`` and
``counterpropagating_weighted_outgoing_total`` must remain zero, while
``counterpropagating_direct_incident_total`` and
``counterpropagating_direct_outgoing_total`` must be positive. The matched-flux
candidate must report ``counterpropagating_matched_no_exchange: True`` in that
case, while the direct-channel candidate must report
``counterpropagating_direct_matched_absorbs_hidden_flux: True`` with negative
exchange power. This compares two candidate laws on the same hidden-flux
configuration; it still does not validate either as a finished absorber.
