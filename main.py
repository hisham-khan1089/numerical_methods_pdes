from pdes import Heat1D, Heat2D, Poisson2D
import numpy as np

heat_1D_solver = Heat1D(length=20, k=1, delta_x=0.5, delta_t=0.1, iterations=500, u_0=100, u_L=100)
heat_1D_solver.animate(filename="heat_1D_sim.mp4")

heat_2D_solver = Heat2D(length=20, k=0.5, delta_x=0.5, delta_t=0.1, iterations=500, u_bottom=100, u_top=100, u_left=100, u_right=100)
heat_2D_solver.animate(filename="heat_2D_sim.mp4")

def source_f(x, y):
    return np.sin(5 * np.pi * x / 100) + np.cos(5 * np.pi * y / 100)

poisson_2D_solver = Poisson2D(length=100, delta_x=1.0, source_func=source_f)
poisson_2D_solver.plot(filename="2D_poisson.png")