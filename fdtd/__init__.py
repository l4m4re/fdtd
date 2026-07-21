""" Python 3D FDTD Simulator """

__author__ = "Floris laporte"
__version__ = "0.3.6"

from .grid import Grid
from .aethergrid import AetherGrid
from .potentialgrid import PotentialGrid, C_LIGHT, VISCOSITY
from .sources import (
    PointSource,
    LineSource,
    PlaneSource,
    AetherPointSource,
    AetherLineSource,
    AetherNativeAngularPointSource,
)
from .detectors import (
    LineDetector,
    BlockDetector,
    CurrentDetector,
    AetherNativeAngularDetector,
)
from .objects import Object, AbsorbingObject, AnisotropicObject
from .boundaries import (
    PeriodicBoundary,
    PML,
    AetherAngularSpongeBoundary,
    AetherAngularNoExchangeBoundary,
    AetherAngularReflectiveBoundary,
    AetherAngularMatchedFluxBoundary,
    AetherAngularDirectMatchedFluxBoundary,
)
from .backend import backend
from .backend import set_backend
from .fourier import FrequencyRoutines
from .visualization import dB_map_2D, plot_detection
