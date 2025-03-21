
from operators import delta_C, delta_I, grad, div, curl_point, curl_surface

'''

Let's implement a new grid class for simuating the fields defined in de delta_C
and delta_L functions in operators.py

The goal is to have a grid that is able to simulate the fields in the delta_C
and delta_L functions in operators.py in a manner similar to the Grid class in
grid.py and LGGrid class in lgrid.py.

The goal of the project is to refactor the existing code from the existing
electromagnetic simulation code into a new simulator for the vector laplacian
derived fields in the delta_C and delta_L functions in operators.py.

The fields are to be evaluated using a staggered grid.

'''


class LGrid:
    def __init__(self, size: Tuple[int, int, int], dt: float):
        """
        Initialize the LGrid with the given size and time step.
        """
        self._size = size
        self._dt = dt
        self._C = bd.zeros((size[0], size[1], size[2], 3))
        self._I = bd.zeros((size[0], size[1], size[2], 3))

    @property
    def C(self) -> Tensorlike:
        """
        Get the C field.
        """
        return self._C

    @property
    def I(self) -> Tensorlike:
        """
        Get the I field.
        """
        return self._I

    def step(self):
        """
        Perform a simulation step.
        """
        self._C += self._dt * delta_I(self._I)
        self._I += self._dt * delta_C(self._C)

    def run(self, steps: int):
        """
        Run the simulation for the given number of steps.
        """
        for _ in range(steps):
            self.step()
            
            
            