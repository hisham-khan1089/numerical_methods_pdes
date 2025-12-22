import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time

# solving the nonhomogeneous case of the 2D heat equation, involving the addition
# of a function f(x, y, t) to the homogeneous case

start = time.perf_counter()

L = 50 # side length of square plate
iterations = 500

k = 1.0
delta_x = 1

delta_t = 0.25
s = (k * delta_t) / (delta_x ** 2)

# generate square plate for each iteration
u = np.empty((iterations, L, L))
u.fill(0)

def f(x, y, t):
    return 50*(np.sin(5*np.pi*x/L))**2 * (np.sin(5*np.pi*y/L))**2 * np.e**(-t)

for i in range(L): # note that range() should be taking L/delta_x, assuming that L/delta_x is an integer
    for j in range(L):
        u[0][j][i] = f(i*delta_x, j*delta_x, 0)

# boundary conditions
u_top = 0.0
u_left = 0.0
u_bottom = 0.0
u_right = 0.0

# # set initial conditions everywhere onto grid
# u_init = 0
# u.fill(u_init)

# set boundary conditions for all iterations
u[:, (L-1):, :] = u_top
u[:, :, :1] = u_left
u[:, :1, 1:] = u_bottom
u[:, :, (L-1):] = u_right

def calculate(u):
    """
    Use partial difference equation for 2D heat equation
    to determine solution for every time iteration in u
    """
    for k in range(0, iterations-1):
        for i in range(1, L-1, delta_x):
            for j in range(1, L-1, delta_x):
                u[k+1][i][j] = s * (u[k][i+1][j] + u[k][i-1][j] + u[k][i][j+1] + u[k][i][j-1] - 4*u[k][i][j]) + u[k][i][j] + \
                delta_t*f(i*delta_x, j*delta_x, k*delta_t)

    return u

x = np.arange(L, step = delta_x)
y = np.arange(L, step = delta_x)
X, Y = np.meshgrid(x, y)

u = calculate(u)

fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))

def update_plot(k, u_k):
    """
    function that generates the plot for a given u_k

    k: frame
    u_k: the temperature profile at frame k
    """
    plt.cla()
    ax.plot_surface(X, Y, u_k, cmap=plt.cm.viridis)

    ax.set_title(f"Temperature of 2D plate at t = {k * delta_t:.3f} s")
    ax.set_xlabel("x position")
    ax.set_ylabel("y position")
    ax.set_zlabel("Temperature (C)")

    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_zlim(0, 100)

    return ax
    
def animate(k):
    """
    function that is passed to FuncAnimation that does the updating
    k: frame
    """
    return update_plot(k, u[k])

animation = FuncAnimation(fig=fig, func=animate, frames = iterations, interval = 40, repeat = False)
animation.save(f'nonhomogeneous_heat_sim_2D.gif')
# plt.show() # NOTE: Execution time will only be printed if you close the figure window 

end = time.perf_counter()
duration = end-start
print(f"Execution time: {duration:.4f} seconds")


