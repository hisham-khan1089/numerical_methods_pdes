import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
import time

class Poisson2D:
    length: float         # square domain side length
    delta_x: float        # spatial discretization
    source_func: callable # source function f(x, y)

    def __init__(self, length: float, delta_x: float, source_func):
        if length <= 0:
            raise ValueError("length must be a positive number.")
        if delta_x <= 0:
            raise ValueError("Spatial discretization (delta_x) must be a positive number.")
        if delta_x >= length:
            raise ValueError("delta_x must be smaller than length.")

        self.length = length
        self.delta_x = delta_x
        self.source_func = source_func
        self.solved = False

    def _initialize_grid(self):
        """Initialize grid dimensions and spacing. (private method)"""
        print("Initializing grid...")
        self.n_x = int(self.length // self.delta_x)
        self.n_y = self.n_x
        self.dx = self.length / self.n_x
        self.dy = self.dx

    def solve(self):
        """Solve the 2D Poisson equation using finite differences and sparse matrix solvers. (public method)"""
        print("Calculating 2D Poisson solution...")
        self.solved = True
        start_time = time.perf_counter()

        self._initialize_grid()

        n_x = self.n_x
        n_y = self.n_y
        dx = self.dx
        dy = self.dy

        A = sparse.coo_array(np.eye(n_x - 1))
        diagonals = [1 / (dy ** 2), -2 / (dx ** 2) - 2 / (dy ** 2), 1 / (dy ** 2)]
        sub_matrix = sparse.diags_array(diagonals, offsets=[-1, 0, 1], shape=(n_y - 1, n_y - 1))

        main_matrix = sparse.kron(A, sub_matrix, format='csr')

        remaining_diags = sparse.diags_array([1 / (dx ** 2), 1 / (dx ** 2)], offsets=[1 - n_y, n_y - 1],
                                             shape=((n_y - 1) * (n_x - 1), (n_y - 1) * (n_x - 1)), format='csr')
        main_matrix = main_matrix + remaining_diags

        f_vals = np.array([self.source_func(i * dx, j * dy) for i in range(1, n_x) for j in range(1, n_y)])
        u_vec = sparse.linalg.spsolve(main_matrix, f_vals)

        self.u = np.reshape(u_vec, (n_x - 1, n_x - 1))

        end_time = time.perf_counter()
        self.solve_time = end_time - start_time
        print(f"Calculation completed! Execution time: {self.solve_time:.4f} s")

    def _generate_and_save_plot(self, filename: str):
        """Generate and save 3D surface plot of the numerical solution. (private method)"""
        print("Plotting results...")
        start_time = time.perf_counter()

        x = np.arange(1, self.length, step=self.dx)
        y = np.arange(1, self.length, step=self.dy)
        X, Y = np.meshgrid(x, y)

        fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
        ax.plot_surface(X, Y, self.u, cmap=plt.cm.viridis)

        ax.set_title("Numerical Solution to 2D Poisson Equation")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("u")

        ax.set_xlim(0, self.length)
        ax.set_ylim(0, self.length)
        ax.view_init(25, -45, 0)

        plt.savefig(filename, dpi=300)

        end_time = time.perf_counter()
        self.plot_time = end_time - start_time
        print(f"Saved plot! Execution time: {self.plot_time:.4f} s")

    def plot(self, filename: str = "2D_poisson.png"):
        """Solve the equation if needed and generate/save the 3D surface plot. (public method)"""
        if not self.solved:
            self.solve()
        self._generate_and_save_plot(filename)


if __name__ == "__main__":
    def source_f(x, y):
        return np.sin(5 * np.pi * x / 100) + np.cos(5 * np.pi * y / 100)

    solver = Poisson2D(length=100, delta_x=1.0, source_func=source_f)
    solver.plot(filename="2D_poisson.png")