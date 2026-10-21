import matplotlib.pyplot as plt
from scipy import sparse
import numpy as np

# dimensions of 2D domain
Lx = 100
Ly = Lx 

# last grid point indices for x- and y-axis (index starts at 0)
n_x = 100
n_y = 100
dx = Lx/n_x
dy = Ly/n_y

def f(x, y):
    """function on right-hand side of Poisson equation (source function)"""
    return np.sin(5*np.pi*x/Lx) + np.cos(5*np.pi*y/Ly)

A = sparse.coo_array(np.eye(n_x-1))
diagonals = [1/(dy**2), -2/(dx**2)-2/(dy**2), 1/(dy**2)]
sub_matrix = sparse.diags_array(diagonals, offsets = [-1,0,1], shape = (n_y-1, n_y-1)) # diagonal submatrix

Matrix = sparse.kron(A, sub_matrix, format = 'csr') # Kronecker tensor product

remaining_diags = sparse.diags_array([1/(dx**2), 1/(dx**2)], offsets=[1-n_y, n_y-1], 
                                     shape = ((n_y-1)*(n_x-1), (n_y-1)*(n_x-1)), format = 'csr')
Matrix = Matrix + remaining_diags

f = np.array([f(i*dx,j*dy) for i in range(1, n_x) for j in range(1, n_y)]) # discretized source function values
u = sparse.linalg.spsolve(Matrix, f) # solve linear system for potential function u

# convert solution from vector to 2D array (so that we can plot it)
u_2D=[]
for i in reversed(range(0, n_x-1)):
    u_2D.append(u[i*(n_x-1):(i+1)*(n_x-1)])
u_2D = np.array(u_2D)

x = np.arange(1, Lx, step = dx)
y = np.arange(1, Ly, step = dy)
X, Y = np.meshgrid(x,y)

fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))

ax.plot_surface(X, Y, u_2D, cmap=plt.cm.viridis)

ax.set_title(f"Numerical Solution to 2D Poisson Equation")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("u")

ax.set_xlim(0, Lx)
ax.set_ylim(0, Ly)
ax.view_init(25, -45, 0)
plt.savefig("2D_poisson.png", dpi=300)
plt.show()

class Poisson2D:
    length: float       # rod length
    k: float            # thermal diffusivity
    delta_x: float      # spatial discretization
    iterations: int     # number of time iterations
    func: function      # function to define initial conditions

    def __init__(self, length, k, delta_x, iterations, u_0, u_L):

        if not isinstance(iterations, int) and iterations > 1:
            raise TypeError("iterations must be a positive integer greater than 1.")
        if length <= 0:
            raise ValueError("length must be a positive number.")
        if delta_x <= 0:
            raise ValueError("Spatial discretizations (delta_x) must be a positive number.")

        self.length = length
        self.k = k
        self.delta_x = delta_x
        self.iterations = iterations
        self.u_0 = u_0
        self.u_L = u_L 

        self.solved = False

    def _initialize_grid(self):
        
