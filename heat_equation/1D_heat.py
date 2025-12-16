import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import time

start = time.perf_counter()

L = 50 # 1D rod length
iterations = 1000

k = 1.0         # thermal diffusivity
delta_x = 1     # space discretization
delta_t = 0.25    # time interval
s = (k * delta_t) / (delta_x ** 2)

if s > 0.5:
    print("WARNING: Stability conditions have been violated")

# boundary conditions
u_0 = 100
u_L = 100

u=np.zeros((iterations, L+1))
u[:, :1] = u_0
u[:, L:] = u_L

def update_temps(u): 
    for j in range(0, iterations-1):
        for i in range(1, len(u[0])-1): 
            u[j+1][i] = s * (u[j][i+1] - 2*u[j][i] + u[j][i-1]) + u[j][i] # partial difference equation
    return u

def plot_temp(u_t, t):
    plt.clf() # clear entire figure

    plt.title(f"Temperature at t = {t * delta_t:.3f} unit time")
    plt.xlabel("x-position")
    plt.ylabel("Temperature")

    plt.plot(u_t)
    plt.xlim(0, L)
    plt.ylim(bottom = 0)

    return plt

def animate(t):
    plot_temp(u[t],t)

u = update_temps(u)
anim = FuncAnimation(plt.figure(), animate, interval=15, frames=iterations, repeat=False)
anim.save(f'heat_sim_1D.mp4')

end = time.perf_counter()
duration = end-start
print(f"Execution time: {duration:.4f} seconds")






