"""Discrete operators used by the educational ``fdtd`` package.

The classic Maxwell solver uses a Yee-style staggering:

- ``E``-type fields live on edge-like slots;
- ``H``-type fields live on face-like slots;
- ``curl_E`` maps edge data to face data; and
- ``curl_H`` maps face data to edge data.

The experimental aether path currently reuses the same discrete patterns with
more neutral names:

- ``curl_edge_to_face`` for circulation from edge-like slots to face-like
  slots;
- ``curl_face_to_edge`` for the reverse map;
- ``linear_to_angular_bridge`` for the current temporary map from the linear
  translational sector into the angular sector;
- ``angular_to_linear_bridge`` for the current temporary map from the angular
  sector back into the linear sector;
- ``apply_angular_metric`` for the meter-carrying length factor that the
  angular sector is expected to need once its native geometry is explicit;
- ``angular_clock_eigenvalue`` for the native two-clock angular benchmark
  before it is coupled into the grid state;
- ``grad`` for cell-scalar to oriented vector differences; and
- ``div`` for contracting vector differences back to a scalar slot.

This file intentionally favors clarity over aggressive abstraction because the
package is meant to be inspectable and educational.
"""

## Imports

# typing
from .typing_ import Tensorlike

# relative
from .backend import backend as bd

## Functions
def curl_E(E: Tensorlike) -> Tensorlike:
    """Map an edge-like field to a face-like field by taking its curl.

    Args:
        E: Electric field sampled on the edge-like slots of the Yee grid.

    Returns:
        The discrete curl of ``E`` on face-like slots.
    """
    curl = bd.zeros(E.shape, dtype=E.dtype)

    curl[:, :-1, :, 0] += E[:, 1:, :, 2] - E[:, :-1, :, 2]
    curl[:, :, :-1, 0] -= E[:, :, 1:, 1] - E[:, :, :-1, 1]

    curl[:, :, :-1, 1] += E[:, :, 1:, 0] - E[:, :, :-1, 0]
    curl[:-1, :, :, 1] -= E[1:, :, :, 2] - E[:-1, :, :, 2]

    curl[:-1, :, :, 2] += E[1:, :, :, 1] - E[:-1, :, :, 1]
    curl[:, :-1, :, 2] -= E[:, 1:, :, 0] - E[:, :-1, :, 0]

    return curl


def curl_edge_to_face(v: Tensorlike) -> Tensorlike:
    """Aether-oriented alias for ``curl_E``.

    The experimental STPT implementation currently interprets a general vector
    field ``v`` as living on edge-like slots and reuses the same discrete
    circulation pattern as the Maxwell ``E -> H`` curl.
    """
    return curl_E(v)



def curl_H(H: Tensorlike) -> Tensorlike:
    """Map a face-like field to an edge-like field by taking its curl.

    Args:
        H: Magnetic field sampled on face-like Yee slots.

    Returns:
        The discrete curl of ``H`` on edge-like slots.
    """
    curl = bd.zeros(H.shape, dtype=H.dtype)

    curl[:, 1:, :, 0] += H[:, 1:, :, 2] - H[:, :-1, :, 2]
    curl[:, :, 1:, 0] -= H[:, :, 1:, 1] - H[:, :, :-1, 1]

    curl[:, :, 1:, 1] += H[:, :, 1:, 0] - H[:, :, :-1, 0]
    curl[1:, :, :, 1] -= H[1:, :, :, 2] - H[:-1, :, :, 2]

    curl[1:, :, :, 2] += H[1:, :, :, 1] - H[:-1, :, :, 1]
    curl[:, 1:, :, 2] -= H[:, 1:, :, 0] - H[:, :-1, :, 0]

    return curl


def curl_face_to_edge(A: Tensorlike) -> Tensorlike:
    """Aether-oriented alias for ``curl_H``.

    This is the reverse circulation map used when a face-like field such as a
    torque or vector-potential-like quantity needs to be expressed back on
    edge-like slots.
    """
    return curl_H(A)


def apply_angular_metric(
    field: Tensorlike,
    metric_length: Tensorlike = None,
) -> Tensorlike:
    """Apply an angular-sector metric length to a field.

    Args:
        field: Angular-sector field with shape ``(Nx, Ny, Nz, 3)``.
        metric_length: Optional metric or lever-arm length with shape
            ``(Nx, Ny, Nz, 1)`` or ``(Nx, Ny, Nz, 3)``.

    Returns:
        ``field`` unchanged when ``metric_length`` is ``None``; otherwise the
        field weighted by the supplied metric length.

    Notes:
        This helper is intentionally simple. It marks where the angular sector
        is expected to carry an extra meter-valued geometric factor without
        prematurely fixing the final constitutive law.
    """

    if metric_length is None:
        return field
    return field * metric_length


def angular_clock_eigenvalue(
    omega_t: Tensorlike,
    omega_p: Tensorlike,
    gamma: Tensorlike,
    delta: float,
) -> Tensorlike:
    """Return the finite-step native angular-clock eigenvalue benchmark.

    Args:
        omega_t: Toroidal or local "around" angular rate ``[1/s]``.
        omega_p: Poloidal or local "through" angular rate ``[1/s]``.
        gamma: Hyperbolic envelope, boost, or dilation rate ``[1/s]``.
        delta: Step size used in the angular-clock finite difference.

    Returns:
        The benchmark eigenvalue
        ``-4/delta**2*sin(delta*omega_t/2)**2
        -4/delta**2*sin(delta*omega_p/2)**2
        +4/delta**2*sinh(delta*gamma/2)**2``.

    Notes:
        This is a pure analytic benchmark for the future native angular sector.
        It deliberately does not update ``AetherGrid`` state or define a charge
        observable.
    """

    if delta == 0:
        raise ValueError("delta must be non-zero")

    omega_t = bd.asarray(omega_t)
    omega_p = bd.asarray(omega_p)
    gamma = bd.asarray(gamma)
    half_delta = 0.5 * delta

    sin_t = bd.sin(half_delta * omega_t)
    sin_p = bd.sin(half_delta * omega_p)
    sinh_h = 0.5 * (
        bd.exp(half_delta * gamma) - bd.exp(-half_delta * gamma)
    )

    scale = 4.0 / (delta * delta)
    return scale * (-sin_t * sin_t - sin_p * sin_p + sinh_h * sinh_h)


def linear_to_angular_bridge(
    linear_field: Tensorlike,
    metric_length: Tensorlike = None,
) -> Tensorlike:
    """Map the linear sector into the current temporary angular bridge.

    Args:
        linear_field: Linear translational state, currently sampled in the same
            array layout as the educational aether prototype.
        metric_length: Optional angular metric length used as a placeholder for
            the meter-valued geometry expected in the native angular sector.

    Returns:
        The current bridge field on angular slots.

    Notes:
        This is still a bridge, not the final STPT angular update law. Without
        ``metric_length`` it reproduces the legacy ``curl_edge_to_face`` path.
    """

    angular_field = curl_edge_to_face(linear_field)
    return apply_angular_metric(angular_field, metric_length)


def angular_to_linear_bridge(
    angular_field: Tensorlike,
    metric_length: Tensorlike = None,
) -> Tensorlike:
    """Map the angular sector back into the current temporary linear bridge.

    Args:
        angular_field: Angular-sector field on face-like slots.
        metric_length: Optional angular metric length used as a simple weighting
            factor in the current bridge implementation.

    Returns:
        The reverse bridge field on linear slots.

    Notes:
        This currently reuses ``curl_face_to_edge`` so the runnable prototype
        stays compatible with the existing baseline while the native angular
        geometry is still being derived.
    """

    weighted_angular_field = apply_angular_metric(angular_field, metric_length)
    return curl_face_to_edge(weighted_angular_field)



def div(v: Tensorlike) -> Tensorlike:
    """Contract an edge-like vector field into a scalar difference field.

    Args:
        v: Vector field with shape ``(Nx, Ny, Nz, 3)``.

    Returns:
        Scalar field with shape ``(Nx, Ny, Nz, 1)``.

    Notes:
        The current implementation uses centered finite differences arranged to
        stay compatible with the package's simple collocated aether prototype.
    """
    div_v = bd.zeros((v.shape[0], v.shape[1], v.shape[2], 1), dtype=v.dtype)

    # Compute the x-component of the divergence
    div_v[1:-1, :, :, 0] += (v[1:-1, :, :, 0] - v[:-2, :, :, 0])
    div_v[:-2, :, :, 0]  -= (v[1:-1, :, :, 0] - v[2:, :, :, 0])

    # Compute the y-component of the divergence
    div_v[:, 1:-1, :, 0] += (v[:, 1:-1, :, 1] - v[:, :-2, :, 1])
    div_v[:, :-2, :, 0]  -= (v[:, 1:-1, :, 1] - v[:, 2:, :, 1])

    # Compute the z-component of the divergence
    div_v[:, :, 1:-1, 0] += (v[:, :, 1:-1, 2] - v[:, :, :-2, 2])
    div_v[:, :, :-2, 0]  -= (v[:, :, 1:-1, 2] - v[:, :, 2:, 2])

    return div_v



def grad(p: Tensorlike) -> Tensorlike:
    """Map a scalar field to oriented finite differences along each axis.

    Args:
        p: Scalar field with shape ``(Nx, Ny, Nz, 1)``.

    Returns:
        Vector field with shape ``(Nx, Ny, Nz, 3)`` whose components represent
        forward differences along ``x``, ``y``, and ``z``.
    """
    grad = bd.zeros((p.shape[0], p.shape[1], p.shape[2], 3), dtype=p.dtype)
    
    # Compute the x-component of the gradient
    grad[:-1, :, :, 0] = (p[1:, :, :, 0] - p[:-1, :, :, 0])
    
    # Compute the y-component of the gradient
    grad[:, :-1, :, 1] = (p[:, 1:, :, 0] - p[:, :-1, :, 0])
    
    # Compute the z-component of the gradient
    grad[:, :, :-1, 2] = (p[:, :, 1:, 0] - p[:, :, :-1, 0])
    
    return grad
