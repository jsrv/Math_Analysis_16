# iterative_solvers.py
"""Volume 1A: Iterative Solvers.
<Juan>
<M345GL>
<Date>
"""

import numpy as np
import random
from scipy import linalg as la
from scipy import sparse
import matplotlib as plt
import time


# Helper function
def diag_dom(n, num_entries=None):
    """Generate a strictly diagonally dominant nxn matrix.

    Inputs:
        n (int): the dimension of the system.
        num_entries (int): the number of nonzero values. Defaults to n^(3/2)-n.

    Returns:
        A ((n,n) ndarray): An nxn strictly diagonally dominant matrix.
    """
    if num_entries is None:
        num_entries = int(n**1.5) - n
    A = np.zeros((n,n))
    rows = np.random.choice(np.arange(0,n), size=num_entries)
    cols = np.random.choice(np.arange(0,n), size=num_entries)
    data = np.random.randint(-4, 4, size=num_entries)
    for i in xrange(num_entries):
        A[rows[i], cols[i]] = data[i]
    for i in xrange(n):
        A[i,i] = np.sum(np.abs(A[i])) + 1
    return A


# Problems 1 and 2
def jacobi_method(A, b, tol=1e-8, maxiters=100, plot=False):
    """Calculate the solution to the system Ax = b voa the Jacobi Method.

    Inputs:
        A ((n,n) ndarray): A square matrix.
        b ((n,) ndarray): A vector of length n.
        tol (float, opt): the convergence tolerance.
        maxiters (int, opt): the maximum number of iterations to perform.
        plot (bool, opt): if True, plot the convergence rate of the algorithm.
            (this is for Problem 2).

    Returns:
        x ((n,) ndarray): the solution to system Ax = b.
    """
    x = np.random.random(len(b))
    A_a = np.copy(A)

    D = np.diag(A_a)
    d_inv= []
    for i in xrange(len(b)):
        d_inv.append(1./D[i])

    d_inv = np.diag(d_inv)

    print x_k
    print d_inv

    this = False
    init_tries= maxiters
    tries = maxiters
    error = []

    while this is False and tries > 0:
        x_kmasuno = x_k + d_inv.dot(b - A.dot(x_k))

        if la.norm((x_k - x_kmasuno), ord=np.inf) <= tol:
            this = True
            differnece = (A.dot( x_kmasuno ) - b)
            error.append( la.norm( differnece, ord=np.inf))

        else:
            this = False
            error.append(A.dot(x_kmasuno)-b)
        tries -= 1

    if plot == True:
        domain = np.arange(0, init_tries - tries)
        plt.title("convergence of the Jacobi Method")

        plt.semilogy(domain, error.append)
        plt.xlabel("Number of iterations", fontsize=14)
        plt.ylabel("Absolute value of approximation", fontsize=14)



# Problem 3
def gauss_seidel(A, b, tol=1e-8, maxiters=100, plot=False):
    """Calculate the solution to the system Ax = b via the Gauss-Seidel Method.

    Inputs:
        A ((n,n) ndarray): A square matrix.
        b ((n,) ndarray): A vector of length n.
        tol (float, opt): the convergence tolerance.
        maxiters (int, opt): the maximum number of iterations to perform.
        plot (bool, opt): if True, plot the convergence rate of the algorithm.

    Returns:
        x ((n,) ndarray): the solution to system Ax = b.
    """

    while this is False and tries > 0:
        xi_kmasuno = xk_i + (1/a_ii)(b_i - A[i].T * xk_i)

        if la.norm((x_k - x_kmasuno), inf) <= tol:
            this = True
            differnece = (A.dot( x_kmasuno ) - b)
            error.append( la.norm( differnece, ord=np.inf))

        else:
            this = False
            error.append(A.dot(x_kmasuno)-b)
        tries -= 1

    if plot == True:
        domain = np.arange(0, init_tries - tries)
        plt.title("convergence of the Gauss-Seidel")

        plt.semilogy(domain, error.append)
        plt.xlabel("Number of iterations", fontsize=14)
        plt.ylabel("Absolute value of approximation", fontsize=14)

# Problem 4
def prob4():
    """For a 5000 parameter system, compare the runtimes of the Gauss-Seidel
    method and la.solve(). Print an explanation of why Gauss-Seidel is so much
    faster.
    """
    domain = 2**np.arange(5,12)
    la_solve = []
    gaus = []

    for i in domain:
        print i
        A = diag_dom(i)
        b = np.random.rand(i,1)

        c = time.time()
        queso = la.solve(A,b)
        la_solve.append(c-time.time())

        c = time.time()
        queso = gauss_seidel(A,b)
        gaus.append(c-time.time())

    plt.loglog(domain, la_solve, 'g.-', linewidth=2, markersize=15, label = "lin solve")
    plt.title("Sparsed Solve", fontsize=18)
    plt.loglog(domain, gaus, 'b.-', linewidth=2, markersize=15, label = "Gauss-Seidel")
    plt.title("Dense Solve", fontsize=18)
    plt.legend

    plt.show()



# Problem 5
def sparse_gauss_seidel(A, b, tol=1e-8, maxiters=100):
    """Calculate the solution to the sparse system Ax = b via the Gauss-Seidel
    Method.

    Inputs:
        A ((n,n) csr_matrix): An nxn sparse CSR matrix.
        b ((n,) ndarray): A vector of length n.
        tol (float, opt): the convergence tolerance.
        maxiters (int, opt): the maximum number of iterations to perform.

    Returns:
        x ((n,) ndarray): the solution to system Ax = b.
    """
    rowstart = A.indptr[i]
    rowend = A.indptr[i+1]
    Aix = np.dot(A.data[rowstart:rowend], x[A.indices[rowstart:rowend]])

    while this is False and tries > 0:
        xi_kmasuno = xk_i + (1/a_ii)(b_i - A[i].T * xk_i)

        if la.norm((x_k - x_kmasuno), inf) <= tol:
            this = True
            differnece = (A.dot( x_kmasuno ) - b)
            error.append( la.norm( differnece, ord=np.inf))

        else:
            this = False
            error.append(A.dot(x_kmasuno)-b)
        tries -= 1

# Problem 6
def sparse_sor(A, b, omega, tol=1e-8, maxiters=100):
    """Calculate the solution to the system Ax = b via Successive Over-
    Relaxation.

    Inputs:
        A ((n,n) csr_matrix): An nxn sparse matrix.
        b ((n,) ndarray): A vector of length n.
        omega (float in [0,1]): The relaxation factor.
        tol (float, opt): the convergence tolerance.
        maxiters (int, opt): the maximum number of iterations to perform.

    Returns:
        x ((n,) ndarray): the solution to system Ax = b.
    """
    raise NotImplementedError("Problem 6 Incomplete")


# Problem 7
def finite_difference(n):
    """Return the A and b described in the finite difference problem that
    solves Laplace's equation.
    """
    raise NotImplementedError("Problem 7 Incomplete")


# Problem 8
def compare_omega():
    """Time sparse_sor() with omega = 1, 1.05, 1.1, ..., 1.9, 1.95, tol=1e-2,
    and maxiters = 1000 using the A and b generated by finite_difference()
    with n = 20. Plot the times as a function of omega.
    """
    raise NotImplementedError("Problem 8 Incomplete")


# Problem 9
def hot_plate(n):
    """Use finite_difference() to generate the system Au = b, then solve the
    system using SciPy's sparse system solver, scipy.sparse.linalg.spsolve().
    Visualize the solution using a heatmap using np.meshgrid() and
    plt.pcolormesh() ("seismic" is a good color map in this case).
    """
    raise NotImplementedError("Problem 9 Incomplete")

if __name__ == "__main__":
    prob4()
"""
    a = diag_dom(4)
    b = np.random.rand(4,1)
    jacobi_method(a,b)
"""
