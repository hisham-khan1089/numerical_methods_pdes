import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import time
from numba import njit

@njit
def initialize_grid(iterations, num_grid_points, u_0, u_L):
    u=np.zeros((iterations+1, num_grid_points+1))
    u[:, :1] = u_0
    u[:, num_grid_points:] = u_L
    return u

@njit
def solve_u(u, s, iterations, num_grid_points):
    for j in range(0, iterations):
        for i in range(1, num_grid_points): 
            u[j+1, i] = s * (u[j, i+1] - 2*u[j, i] + u[j, i-1]) + u[j, i]
    return u

class heat_1D:
    length: float       # rod length
    k: float            # thermal diffusivity
    delta_x: float      # spatial discretization
    delta_t: float      # temporal discretization
    iterations: int     # number of time iterations
    u_0: float          # left boundary condition (x=0)
    u_L: float          # right boundary condition (x=L)

    def __init__(self, length, k, delta_x, delta_t, iterations, u_0, u_L):

        if (k * delta_t) / (delta_x ** 2) >= 0.5:
            raise ValueError("Stability condition (k * delta_t) / (delta_x^2) < 0.5 has been violated.")
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
        self.u_0 = u_0
        self.u_L = u_L 

        self.solved = False

    def _initialize_grid(self):
        """Initialize the temperature grid for all time iterations. (private method)"""
        print("Initializing grid...")
        self.num_grid_points = int(self.length // self.delta_x)
        self.u = initialize_grid(self.iterations, self.num_grid_points, self.u_0, self.u_L)

    def solve(self): 
        """Solve the temperature grid for each time iteration. (public method)"""

        print("Calculating temperature grid values...")
        self.solved = True
        start_time = time.perf_counter()

        self._initialize_grid()
        s = (self.k * self.delta_t) / (self.delta_x **2)
        self.u = solve_u(self.u, s, self.iterations, self.num_grid_points)

        end_time = time.perf_counter()
        self.solve_time = end_time - start_time
        print(f"Calculation completed! Execution time: {self.solve_time:.4f} s")

    def _generate_and_save_animation(self, filename: str):
        """Generate mp4 file containing time evolution of temperature of rod over time. (private method)"""

        print("Animating results...")

        start_time = time.perf_counter()

        fig, ax = plt.subplots()
        ax.set_xlabel("x")
        ax.set_ylabel("Temperature")
        ax.set_xlim(0, self.length)

        min_ylim = min(self.u_0, self.u_L, min(self.u[0]))
        max_ylim = max(self.u_0, self.u_L)
        padding = 0.05 * (max_ylim-min_ylim)
        ax.set_ylim(min_ylim, max_ylim+padding)

        line, = ax.plot([self.delta_x*i for i in range(self.num_grid_points+1)], self.u[0])

        title_text = ax.set_title("")

        def plot_temp(u_t, t):
            title_text.set_text(f"Temperature at t = {t * self.delta_t:.3f} s")
            line.set_ydata(u_t)
            return line, title_text

        def _animate(t):
            return plot_temp(self.u[t],t)

        anim = FuncAnimation(fig, _animate, interval=15, frames=self.iterations+1, repeat=False, blit=True)
        anim.save(filename)

        end_time = time.perf_counter()
        self.animation_time = end_time - start_time
        print(f"Saved animation! Execution time: {self.animation_time:.4f} s")

    def animate(self, filename: str):
        """Solve for time evolution of temperature and create animation. (public method)"""
        if self.solved is False:
            self.solve()
        self._generate_and_save_animation(filename)


if __name__ == "__main__":
    solver = heat_1D(length=20, k=1, delta_x=0.5, delta_t=0.1, iterations=500, u_0=100, u_L=100)
    solver.animate(filename="heat_1D_test.mp4")
