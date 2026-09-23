import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import time
from numba import njit, prange

@njit(parallel=True)
def initialize_grid(iterations, num_grid_points, u_top, u_left, u_bottom, u_right):
    u = np.zeros((iterations+1, num_grid_points+1, num_grid_points+1))
    u[:, (num_grid_points):, :] = u_top
    u[:, :, :1] = u_left
    u[:, :1, 1:] = u_bottom
    u[:, :, (num_grid_points):] = u_right
    return u

@njit(parallel=True, fastmath=True)
def solve_u(u, s, iterations, num_grid_points):
    for k in range(0, iterations):
        for i in prange(1, num_grid_points):
            for j in range(1, num_grid_points):
                u[k+1, i, j] = s * (u[k, i+1, j] + u[k, i-1, j] + u[k, i, j+1] + u[k, i, j-1] - 4*u[k, i, j]) + u[k, i, j]
    return u

class heat_2D:
    length: float       # square surface side length
    k: float            # thermal diffusivity
    delta_x: float      # spacial discretization
    delta_t: float      # temporal discretization
    iterations: float   # number of time iterations
    u_top: float        # top boundary condition (y=L)
    u_bottom: float     # bottom boundary condition (y=0)
    u_left: float       # left boundary condition (x=0)
    u_right: float      # right boundary condition (x=L)

    def __init__(self, length, k, delta_x, delta_t, iterations, u_top, u_bottom, u_left, u_right):

        if (k * delta_t) / (delta_x ** 2) >= 0.25:
            raise ValueError("Stability condition (k * delta_t) / (delta_x^2) < 0.25 has been violated.")
        if not isinstance(iterations, int) and iterations > 1:
            raise TypeError("iterations must be a positive integer greater than 1.")
        if length <= 0:
            raise ValueError("length must be a positive number.")
        if delta_t <= 0 or delta_x <= 0:
            raise ValueError("Spatial and temporal discretizations (delta_x and delta_t) must be positive numbers.")

        self.length = length
        self.k = k
        self.delta_x = delta_x
        self.delta_t = delta_t
        self.iterations = iterations
        self.u_top = u_top
        self.u_bottom = u_bottom
        self.u_left = u_left
        self.u_right = u_right

        self.solved = False

    def _initialize_grid(self):
        """Initialize the temperature grid for all time iterations. (private method)"""
        self.num_grid_points = int(self.length // self.delta_x)
        self.length = self.delta_x * self.num_grid_points
        self.u = initialize_grid(self.iterations, self.num_grid_points, 
                                 self.u_top, self.u_left, self.u_bottom, self.u_right)

    def solve(self):

        self.solved = True
        start_time = time.perf_counter()

        self._initialize_grid()
        s = (self.k * self.delta_t) / (self.delta_x ** 2)
        self.u = solve_u(self.u, s, self.iterations, self.num_grid_points)

        end_time = time.perf_counter()
        self.solve_time = end_time - start_time

    def _generate_and_save_animation(self, filename: str):

        fig, ax = plt.subplots()
        ax.set_xlabel("x")
        ax.set_ylabel("y")

        image = ax.imshow(self.u[0], 
                          cmap='inferno', 
                          vmin=np.min(self.u), vmax=np.max(self.u), 
                          origin='lower',
                          extent=[0, self.length, 0, self.length])
        fig.colorbar(image, ax=ax, label="Temperature")

        title_text = ax.set_title("")

        def plot_temp(u_t, t):
            title_text.set_text(f"Temperature at t = {t * self.delta_t:.3f} s")
            image.set_data(u_t)
            return image, title_text

        def _animate(t):
            return plot_temp(self.u[t],t)

        anim = FuncAnimation(fig, _animate, interval=15, frames=self.iterations+1, repeat=False, blit=True)
        anim.save(filename)
        print("Saved animation!")

    def animate(self, filename: str):
        """Solve for time evolution of temperature and create animation. (public method)"""
        if self.solved is False:
            self.solve()
        self._generate_and_save_animation(filename)


if __name__ == "__main__":
    solver = heat_2D(length=20, k=0.5, delta_x=0.5, delta_t=0.1, iterations=500, u_bottom=100, u_top=100, u_left=100, u_right=100)
    solver.animate(filename="heat_2D_test.mp4")