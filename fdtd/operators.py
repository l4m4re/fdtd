"""Discrete operator helpers for the ``fdtd`` package.

The original Maxwell solver in :mod:`fdtd.grid` uses a Yee-style staggering:

- ``E``-type fields are stored on edge-like slots;
- ``H``-type fields are stored on face-like slots;
- ``curl_E`` maps edge-like data to face-like data;
- ``curl_H`` maps face-like data to edge-like data.

The helpers in this module keep those same discrete patterns but provide
slightly more theory-neutral names for the experimental ``AetherGrid`` path.
They should therefore be read as collocation-aware finite differences rather
than as abstract continuum operators.

The current aether implementation is still experimental. Its long-term intent
is to separate linear and angular sectors more clearly, but for now these
operators remain close to the educational Maxwell discretization so that the
scene-construction API and backend behavior stay familiar.
"""

## Imports

# typing
from .typing_ import Tensorlike

# relative
from .backend import backend as bd

## Functions
def curl_E(E: Tensorlike) -> Tensorlike:
    """Transform an edge-like field into a face-like field via a discrete curl.

    Args:
        E: Electric field stored with the same collocation as the standard
            ``Grid.E`` field.

    Returns:
        An ``H``-type field living on the complementary face-like slots.
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
    """Discrete curl from edge-like slots to face-like slots.

    This is the collocation-neutral alias used by :class:`fdtd.aethergrid.AetherGrid`.
    The stencil is identical to :func:`curl_E`.

    Args:
        v: Vector field whose components are stored on edge-like locations.

    Returns:
        A face-like vector field representing the discrete circulation of ``v``.
    """
    return curl_E(v)



def curl_H(H: Tensorlike) -> Tensorlike:
    """Transform a face-like field into an edge-like field via a discrete curl.

    Args:
        H: Magnetic field stored with the same collocation as the standard
            ``Grid.H`` field.

    Returns:
        An ``E``-type field living on the complementary edge-like slots.
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
    """Discrete curl from face-like slots to edge-like slots.

    This is the collocation-neutral alias used by :class:`fdtd.aethergrid.AetherGrid`.
    The stencil is identical to :func:`curl_H`.

    Args:
        A: Vector field whose components are stored on face-like locations.

    Returns:
        An edge-like vector field representing the discrete circulation of ``A``.
    """
    return curl_H(A)



def div(v: Tensorlike) -> Tensorlike:
    """Compute a discrete divergence-like contraction of a vector field.

    The current ``AetherGrid`` stores its vector fields in collocated
    ``(Nx, Ny, Nz, 3)`` arrays even though the longer-term theory likely
    requires a clearer separation between linear and angular slots. This helper
    therefore uses the package's current experimental convention rather than a
    finalized STPT geometry.

    Args:
        v: Vector field with final axis of length 3.

    Returns:
        Scalar field of shape ``(Nx, Ny, Nz, 1)``.
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
    """Compute a discrete gradient-like map from scalar to vector field.

    Args:
        p: Scalar field with trailing singleton axis, typically
            ``(Nx, Ny, Nz, 1)`` in the current experimental aether path.

    Returns:
        Vector field with shape ``(Nx, Ny, Nz, 3)``.
    """
    grad = bd.zeros((p.shape[0], p.shape[1], p.shape[2], 3), dtype=p.dtype)
    
    # Compute the x-component of the gradient
    grad[:-1, :, :, 0] = (p[1:, :, :, 0] - p[:-1, :, :, 0])
    
    # Compute the y-component of the gradient
    grad[:, :-1, :, 1] = (p[:, 1:, :, 0] - p[:, :-1, :, 0])
    
    # Compute the z-component of the gradient
    grad[:, :, :-1, 2] = (p[:, :, 1:, 0] - p[:, :, :-1, 0])
    
    return grad
