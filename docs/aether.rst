Experimental Aether Grid
========================

The fork in this repository adds an experimental ``AetherGrid`` alongside the
original Maxwell ``Grid``.

Why this exists
---------------

The upstream package is already a strong educational simulator:

- scenes are built by registering sources, objects, boundaries, and detectors;
- examples exist and are easy to compare;
- the backend abstraction already supports NumPy, Torch, and CUDA paths.

That makes it a natural place to explore an aether-style extension without
inventing a completely separate simulator culture.

What ``AetherGrid`` is
----------------------

``AetherGrid`` is an opt-in experimental grid. It is not intended to replace
the classic electromagnetic solver. Instead, it provides a theory-specific path
that tries to evolve an aether velocity field and its derived quantities while
keeping the package's public style familiar.

The current field chain is:

.. code-block:: text

   v -> p, omega -> f, tau -> E, H -> A -> a
     -> dpdt, alpha -> dEdt, dHdt -> dAdt -> j

In words:

- ``v`` is the current aether velocity-like state;
- ``p`` is a scalar potential-like quantity derived by divergence;
- ``omega`` is a rotational quantity derived by curl;
- ``f`` and ``tau`` are force-density-like and torque-density-like quantities;
- ``E`` and ``H`` are derived electromagnetic analogues;
- ``A`` is a vector-potential-like field;
- ``a`` and ``j`` are acceleration and jerk.

Current limitations
-------------------

This path is still experimental and should be read with caution:

- most fields are still collocated in ``(Nx, Ny, Nz, 3)`` arrays;
- the code currently reflects an older constant-``k`` formulation;
- the longer-term STPT program likely needs a clearer separation between
  linear and angular sectors at each point in space.

So ``AetherGrid`` is best understood as an integration surface for ongoing
theory work, not yet as a validated physics solver.

How the operator language should be read
----------------------------------------

The original ``fdtd`` package uses Yee-style staggering:

- ``E``-type fields are edge-like;
- ``H``-type fields are face-like;
- ``curl_E`` maps edge-like data to face-like data;
- ``curl_H`` maps face-like data to edge-like data.

The experimental operator names in :mod:`fdtd.operators` are collocation-aware
aliases of those same discrete patterns:

- ``curl_edge_to_face`` uses the same stencil as ``curl_E``;
- ``curl_face_to_edge`` uses the same stencil as ``curl_H``;
- ``grad`` maps a scalar field to a vector field;
- ``div`` contracts a vector field to a scalar field.

For the long-term STPT direction, these should eventually be reinterpreted as
maps between linear and angular sectors rather than as renamed Maxwell curls.

Public-package guidance
-----------------------

If this experimental path is developed further, the safest package direction is:

1. keep the original ``Grid`` behavior unchanged;
2. keep aether features opt-in rather than implicit;
3. document collocation and operator meaning as clearly as the Maxwell path;
4. compare new predictions against existing Maxwell examples whenever possible.

API reference
-------------

.. automodule:: fdtd.aethergrid
   :members:
   :undoc-members:
   :show-inheritance:
