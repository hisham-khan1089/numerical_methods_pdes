import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import time

class heat_2D:
    length: float       # square surface side length
    k: float            # thermal diffusivity
    delta_x: float      # spacial discretization in x-axis
    delta_y: float      # spacial discretization in y-axis
    delta_t: float      # temporal discretization
    iterations: float   # number of time iterations
    u_top: float        # top boundary condition (y=L)
    u_bottom: float     # bottom boundary condition (y=0)

    def __init__(self, length, k, delta_x, delta_y, delta_t, iterations, u_top, u_bottom, u_left, u_right):

        if (k * delta_t) / (delta_x ** 2) > 0.25:
            raise ValueError("Stability condition (k * delta_t) / (delta_x^2) <= 0.25 has been violated.")
        if not isinstance(iterations, int):
            raise TypeError("iterations must be an integer.")

        self.length = length
        self.k = k
        self.delta_x = delta_x
        self.delta_y = delta_y
        self.delta_t = delta_t
        self.iterations = iterations
        self.u_top = u_top
        self.u_bottm = u_bottom
        self.u_left = u_left
        self.u_right = u_right

    def _initialize_grid(self):
        self.num_grid_x = int(np.floor(self.length/self.delta_x))
        self.num_grid_y = int(np.floor())

        u = np.empty((self.iterations, L/delta_x, L/delta_x))

        # initial conditions everywhere on plate
        u_init = 0.0

        # boundary conditions
        u_top = 100.0
        u_left = 100.0
        u_bottom = 100.0
        u_right = 100.0

        # set initial conditions onto grid
        u.fill(u_init)
