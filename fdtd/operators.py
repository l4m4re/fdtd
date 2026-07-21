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

import numpy as np

# typing
from typing import Tuple

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
        gamma: Independent hyperbolic envelope, boost, or dilation rate
            ``[1/s]`` for the current benchmark. It is not derived from
            ``omega_t`` and ``omega_p`` by this helper.
        delta: Step size used in the angular-clock finite difference.

    Returns:
        The benchmark eigenvalue
        ``-4/delta**2*sin(delta*omega_t/2)**2
        -4/delta**2*sin(delta*omega_p/2)**2
        +4/delta**2*sinh(delta*gamma/2)**2``.

    Notes:
        This is a pure analytic benchmark for the future native angular sector.
        It deliberately does not update ``AetherGrid`` state, infer ``gamma``
        from the two angular clocks, or define a charge observable.
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
        The current implementation uses a simple local finite-difference
        stencil arranged to stay compatible with the package's collocated
        aether prototype.
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


## Legacy staggered-grid compatibility operators
#
# ``PotentialGrid`` uses the earlier STPT staggered layout with separate
# component arrays. Keep these helpers here so legacy simulator knowledge can
# live inside the active ``extern/fdtd`` package while the current ``AetherGrid``
# API continues to use the compact 4D operators above.


def gradient(
    scalar_field: Tensorlike,
    dx: float = 1.0,
    dy: float = None,
    dz: float = None,
) -> Tuple[Tensorlike, Tensorlike, Tensorlike]:
    """Calculate a cell-centered scalar gradient on staggered faces."""

    if dy is None:
        dy = dx
    if dz is None:
        dz = dx

    gx = bd.zeros((scalar_field.shape[0] + 1, scalar_field.shape[1], scalar_field.shape[2]))
    gy = bd.zeros((scalar_field.shape[0], scalar_field.shape[1] + 1, scalar_field.shape[2]))
    gz = bd.zeros((scalar_field.shape[0], scalar_field.shape[1], scalar_field.shape[2] + 1))

    gx[1:-1, :, :] = (scalar_field[1:, :, :] - scalar_field[:-1, :, :]) / dx
    gy[:, 1:-1, :] = (scalar_field[:, 1:, :] - scalar_field[:, :-1, :]) / dy
    gz[:, :, 1:-1] = (scalar_field[:, :, 1:] - scalar_field[:, :, :-1]) / dz

    gx[0, :, :] = gx[1, :, :]
    gx[-1, :, :] = gx[-2, :, :]
    gy[:, 0, :] = gy[:, 1, :]
    gy[:, -1, :] = gy[:, -2, :]
    gz[:, :, 0] = gz[:, :, 1]
    gz[:, :, -1] = gz[:, :, -2]

    return gx, gy, gz


def divergence(
    vx: Tensorlike,
    vy: Tensorlike,
    vz: Tensorlike,
    dx: float = 1.0,
) -> Tensorlike:
    """Calculate divergence from separate face-centered vector components."""

    return (
        vx[1:, :, :] - vx[:-1, :, :]
        + vy[:, 1:, :] - vy[:, :-1, :]
        + vz[:, :, 1:] - vz[:, :, :-1]
    ) / dx


def _curl_face_to_edge_staggered(
    vx: Tensorlike,
    vy: Tensorlike,
    vz: Tensorlike,
    dx: float = 1.0,
) -> Tuple[Tensorlike, Tensorlike, Tensorlike]:
    """Curl of separate face-centered components, returned on edge slots."""

    nx = vx.shape[0] - 1
    ny = vy.shape[1] - 1
    nz = vz.shape[2] - 1

    curl_x = bd.zeros((nx, ny + 1, nz + 1))
    curl_y = bd.zeros((nx + 1, ny, nz + 1))
    curl_z = bd.zeros((nx + 1, ny + 1, nz))

    curl_x[:, 1:-1, 1:-1] = (
        (vz[:, 1:, 1:-1] - vz[:, :-1, 1:-1]) / dx
        - (vy[:, 1:-1, 1:] - vy[:, 1:-1, :-1]) / dx
    )
    curl_y[1:-1, :, 1:-1] = (
        (vx[1:-1, :, 1:] - vx[1:-1, :, :-1]) / dx
        - (vz[1:, :, 1:-1] - vz[:-1, :, 1:-1]) / dx
    )
    curl_z[1:-1, 1:-1, :] = (
        (vy[1:, 1:-1, :] - vy[:-1, 1:-1, :]) / dx
        - (vx[1:-1, 1:, :] - vx[1:-1, :-1, :]) / dx
    )

    curl_x[:, 0, :] = curl_x[:, 1, :]
    curl_x[:, -1, :] = curl_x[:, -2, :]
    curl_x[:, :, 0] = curl_x[:, :, 1]
    curl_x[:, :, -1] = curl_x[:, :, -2]
    curl_y[0, :, :] = curl_y[1, :, :]
    curl_y[-1, :, :] = curl_y[-2, :, :]
    curl_y[:, :, 0] = curl_y[:, :, 1]
    curl_y[:, :, -1] = curl_y[:, :, -2]
    curl_z[0, :, :] = curl_z[1, :]
    curl_z[-1, :, :] = curl_z[-2, :]
    curl_z[:, 0, :] = curl_z[1, :]
    curl_z[:, -1, :] = curl_z[-2, :]
    curl_z[:, :, 0] = curl_z[:, :, 1]
    curl_z[:, :, -1] = curl_z[:, :, -2]

    return curl_x, curl_y, curl_z


def _curl_edge_to_face_staggered(
    wx: Tensorlike,
    wy: Tensorlike,
    wz: Tensorlike,
    dx: float = 1.0,
) -> Tuple[Tensorlike, Tensorlike, Tensorlike]:
    """Curl of separate edge-centered components, returned on face slots."""

    nx_wx, ny_wx, nz_wx = wx.shape
    nx_wy, ny_wy, nz_wy = wy.shape
    nx_wz, ny_wz, nz_wz = wz.shape

    curl_x = bd.zeros((nx_wx + 1, ny_wx, nz_wx))
    curl_y = bd.zeros((nx_wy, ny_wy + 1, nz_wy))
    curl_z = bd.zeros((nx_wz, ny_wz, nz_wz + 1))

    for i in range(1, nx_wx):
        for j in range(1, ny_wx - 1):
            for k in range(1, nz_wx - 1):
                curl_x[i, j, k] = (
                    (wz[i - 1, j, k] - wz[i - 1, j - 1, k]) / dx
                    - (wy[i - 1, j - 1, k] - wy[i - 1, j - 1, k - 1]) / dx
                )

    for i in range(1, nx_wy - 1):
        for j in range(1, ny_wy):
            for k in range(1, nz_wy - 1):
                curl_y[i, j, k] = (
                    (wx[i - 1, j - 1, k] - wx[i - 1, j - 1, k - 1]) / dx
                    - (wz[i, j - 1, k] - wz[i - 1, j - 1, k]) / dx
                )

    for i in range(1, nx_wz - 1):
        for j in range(1, ny_wz - 1):
            for k in range(1, nz_wz):
                curl_z[i, j, k] = (
                    (wy[i, j - 1, k - 1] - wy[i - 1, j - 1, k - 1]) / dx
                    - (wx[i - 1, j, k - 1] - wx[i - 1, j - 1, k - 1]) / dx
                )

    curl_x[0, :, :] = curl_x[1, :, :]
    curl_x[-1, :, :] = curl_x[-2, :, :]
    curl_x[:, 0, :] = curl_x[:, 1, :]
    curl_x[:, -1, :] = curl_x[:, -2, :]
    curl_x[:, :, 0] = curl_x[:, :, 1]
    curl_x[:, :, -1] = curl_x[:, :, -2]
    curl_y[0, :, :] = curl_y[1, :, :]
    curl_y[-1, :, :] = curl_y[-2, :, :]
    curl_y[:, 0, :] = curl_y[:, 1, :]
    curl_y[:, -1, :] = curl_y[:, -2, :]
    curl_y[:, :, 0] = curl_y[:, :, 1]
    curl_y[:, :, -1] = curl_y[:, :, -2]
    curl_z[0, :, :] = curl_z[1, :]
    curl_z[-1, :, :] = curl_z[-2, :]
    curl_z[:, 0, :] = curl_z[1, :]
    curl_z[:, -1, :] = curl_z[-2, :]
    curl_z[:, :, 0] = curl_z[:, :, 1]
    curl_z[:, :, -1] = curl_z[:, :, -2]

    return curl_x, curl_y, curl_z


def curl_face_to_edge(*args):
    """Dispatch curl from face slots to edge slots for current or legacy APIs."""

    if len(args) == 1:
        return curl_H(args[0])
    if len(args) in (3, 4):
        dx = args[3] if len(args) == 4 else 1.0
        return _curl_face_to_edge_staggered(args[0], args[1], args[2], dx)
    raise TypeError("curl_face_to_edge expects 1 argument or vx, vy, vz[, dx]")


def curl_edge_to_face(*args):
    """Dispatch curl from edge slots to face slots for current or legacy APIs."""

    if len(args) == 1:
        return curl_E(args[0])
    if len(args) in (3, 4):
        dx = args[3] if len(args) == 4 else 1.0
        return _curl_edge_to_face_staggered(args[0], args[1], args[2], dx)
    raise TypeError("curl_edge_to_face expects 1 argument or wx, wy, wz[, dx]")


def vector_laplacian(
    vx: Tensorlike,
    vy: Tensorlike,
    vz: Tensorlike,
    dx: float = 1.0,
) -> Tuple[Tensorlike, Tensorlike, Tensorlike]:
    """Calculate the vector Laplacian for separate staggered components."""

    lap_vx = bd.zeros_like(vx)
    lap_vy = bd.zeros_like(vy)
    lap_vz = bd.zeros_like(vz)

    lap_vx[1:-1, 1:-1, 1:-1] = (
        vx[2:, 1:-1, 1:-1]
        + vx[:-2, 1:-1, 1:-1]
        + vx[1:-1, 2:, 1:-1]
        + vx[1:-1, :-2, 1:-1]
        + vx[1:-1, 1:-1, 2:]
        + vx[1:-1, 1:-1, :-2]
        - 6 * vx[1:-1, 1:-1, 1:-1]
    ) / (dx * dx)
    lap_vy[1:-1, 1:-1, 1:-1] = (
        vy[2:, 1:-1, 1:-1]
        + vy[:-2, 1:-1, 1:-1]
        + vy[1:-1, 2:, 1:-1]
        + vy[1:-1, :-2, 1:-1]
        + vy[1:-1, 1:-1, 2:]
        + vy[1:-1, 1:-1, :-2]
        - 6 * vy[1:-1, 1:-1, 1:-1]
    ) / (dx * dx)
    lap_vz[1:-1, 1:-1, 1:-1] = (
        vz[2:, 1:-1, 1:-1]
        + vz[:-2, 1:-1, 1:-1]
        + vz[1:-1, 2:, 1:-1]
        + vz[1:-1, :-2, 1:-1]
        + vz[1:-1, 1:-1, 2:]
        + vz[1:-1, 1:-1, :-2]
        - 6 * vz[1:-1, 1:-1, 1:-1]
    ) / (dx * dx)

    for arr in (lap_vx, lap_vy, lap_vz):
        arr[0, :, :] = arr[1, :, :]
        arr[-1, :, :] = arr[-2, :, :]
        arr[:, 0, :] = arr[:, 1, :]
        arr[:, -1, :] = arr[:, -2, :]
        arr[:, :, 0] = arr[:, :, 1]
        arr[:, :, -1] = arr[:, :, -2]

    return lap_vx, lap_vy, lap_vz


def scalar_laplacian(
    scalar_field: Tensorlike,
    dx: float = 1.0,
    dy: float = None,
    dz: float = None,
) -> Tensorlike:
    """Calculate a scalar Laplacian on cell-centered data."""

    if dy is None:
        dy = dx
    if dz is None:
        dz = dx

    lap = bd.zeros_like(scalar_field)
    lap[1:-1, 1:-1, 1:-1] = (
        (scalar_field[2:, 1:-1, 1:-1] + scalar_field[:-2, 1:-1, 1:-1]) / (dx * dx)
        + (scalar_field[1:-1, 2:, 1:-1] + scalar_field[1:-1, :-2, 1:-1]) / (dy * dy)
        + (scalar_field[1:-1, 1:-1, 2:] + scalar_field[1:-1, 1:-1, :-2]) / (dz * dz)
        - 2
        * scalar_field[1:-1, 1:-1, 1:-1]
        * (1 / (dx * dx) + 1 / (dy * dy) + 1 / (dz * dz))
    )

    lap[0, :, :] = lap[1, :, :]
    lap[-1, :, :] = lap[-2, :, :]
    lap[:, 0, :] = lap[:, 1, :]
    lap[:, -1, :] = lap[:, -2, :]
    lap[:, :, 0] = lap[:, :, 1]
    lap[:, :, -1] = lap[:, :, -2]

    return lap


def helmholtz_decomposition(
    vx: Tensorlike,
    vy: Tensorlike,
    vz: Tensorlike,
    dx: float = 1.0,
) -> Tuple[Tensorlike, Tensorlike, Tensorlike, Tensorlike, Tensorlike, Tensorlike]:
    """Legacy approximate Helmholtz decomposition for separate components."""

    div_v = divergence(vx, vy, vz, dx)
    scalar_potential = bd.zeros_like(div_v)

    for _ in range(30):
        lap_phi = scalar_laplacian(scalar_potential, dx)
        scalar_potential += 0.1 * (div_v - lap_phi)

    vx_irrot, vy_irrot, vz_irrot = gradient(scalar_potential, dx)
    vx_sol = vx - vx_irrot
    vy_sol = vy - vy_irrot
    vz_sol = vz - vz_irrot

    return vx_irrot, vy_irrot, vz_irrot, vx_sol, vy_sol, vz_sol


def quantize_circulation(
    field: Tensorlike,
    dx: float,
    kinematic_viscosity: float,
) -> Tensorlike:
    """Quantize a field component to multiples of a circulation quantum."""

    del dx
    quantum = kinematic_viscosity
    return bd.round(field / quantum) * quantum
