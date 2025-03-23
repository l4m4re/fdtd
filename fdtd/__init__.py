""" Python 3D FDTD Simulator with Space Time Potential Theory extensions """

__author__ = "Floris laporte"
__version__ = "0.3.6"

# Original FDTD components
from .grid import Grid
from .sources import PointSource, LineSource, PlaneSource
from .detectors import LineDetector, BlockDetector, CurrentDetector
from .objects import Object, AbsorbingObject, AnisotropicObject
from .boundaries import PeriodicBoundary, PML
from .backend import backend
from .backend import set_backend
from .fourier import FrequencyRoutines
from .visualization import dB_map_2D, plot_detection

# Space Time Potential Theory extensions
from .ai_aethergrid import AetherGrid, VISCOSITY, C_LIGHT, BACKGROUND_DENSITY
from .ai_operators import (
    grad, div, curl_face_to_edge, curl_edge_to_face, vector_laplacian,
    quantize_circulation
)
