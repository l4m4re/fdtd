class LaplaceGrid:
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
        p = div(self._C)
        L = grad(p)
        A = curl_point(self._C)
        R = curl_surface(A)
        self._C += self._dt * (L - R)

        t = div(self._I)
        Y_l = grad(t)
        W = curl_point(self._I)
        Y_a = curl_surface(W)
        self._I += self._dt * (Y_l - Y_a)

    def run(self, steps: int):
        """
        Run the simulation for the given number of steps.
        """
        for _ in range(steps):
            self.step()