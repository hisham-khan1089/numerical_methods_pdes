import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import time

class heat_1D:
    length: float       # rod length
    k: float            # thermal diffusivity
    delta_x: float      # spacial discretization
    delta_t: float      # temporal discretization
    iterations: int   # number of time iterations
    u_0: float          # left boundary condition (x=0)
    u_L: float          # right boundary condition (x=L)

    def __init__(self, length, k, delta_x, delta_t, iterations, u_0, u_L):

        if (k * delta_t) / (delta_x ** 2) > 0.5:
            raise ValueError("Stability condition (k * delta_t) / (delta_x^2) <= 0.5 has been violated.")
        if not isinstance(iterations, int):
            raise TypeError("iterations must be an integer")

        self.length = length
        self.k = k
        self.delta_x = delta_x
        self.delta_t = delta_t
        self.iterations = iterations
        self.u_0 = u_0
        self.u_L = u_L 

    def __initialize_grid(self):
        """Initialize the temperature grid for all time iterations. (private method)"""
        self.num_grid_points = np.floor(self.length / self.delta_x)
        self.num_grid_points = int(self.num_grid_points)
        self.u=np.zeros((self.iterations, self.num_grid_points+1))
        self.u[:, :1] = self.u_0
        self.u[:, self.num_grid_points:] = self.u_L

    def __solve(self): 
        """Solve the temperature grid for each time iteration. (private method)"""
        self.__initialize_grid()
        s = (self.k * self.delta_t) / (self.delta_x **2)
        print("Solving partial difference equation...")
        for j in range(0, self.iterations-1):
            for i in range(1, len(self.u[0])-1): 
                self.u[j+1, i] = s * (self.u[j, i+1] - 2*self.u[j, i] + self.u[j, i-1]) + self.u[j, i] # partial difference equation
        print("Solving complete!")
        return self

    def __generate_and_save_animation(self, path: str= None):
        """Generate mp4 file containing time evolution of temperature of rod over time. (private method)"""
        print("Animating results...")
        def plot_temp(u_t, t):
            plt.clf() # clear entire figure

            plt.title(f"Temperature at t = {t * self.delta_t:.3f} unit time")
            plt.xlabel("x-position")
            plt.ylabel("Temperature")

            plt.plot([self.delta_x*i for i in range(self.num_grid_points+1)], u_t)
            plt.xlim(0, self.length)
            plt.ylim(bottom = 0)

            return plt

        def animate(t):
            plot_temp(self.u[t],t)

        anim = FuncAnimation(plt.figure(), animate, interval=15, frames=self.iterations, repeat=False)
        anim.save(f'heat_sim_1D.mp4')
        print("Saved results!")

    def solve_and_generate(self):
        """Solve for time evolution of temperature and create animation. (public method)"""
        self.__solve()
        self.__generate_and_save_animation()


