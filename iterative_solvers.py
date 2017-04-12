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
from matplotlib import pyplot as plt
import time
from scipy.sparse import linalg as spla

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
    x_k = np.zeros(len(b))
    A_a = np.copy(A)
    D = np.diag(A_a)
    d_inv= []
    for i in xrange(len(b)):
        d_inv.append(1./D[i])

    d_inv = np.diag(d_inv)
    this = False
    tries = maxiters
    error = []

    while this is False and tries > 0:
        x_kmasuno = x_k + d_inv.dot(b - A.dot(x_k))

        if ((la.norm((x_k - x_kmasuno), ord=np.inf)) <= tol):
            this = True
            difference = (A.dot( x_kmasuno ) - b)
            error.append( la.norm( difference, ord=np.inf))


        else:
            difference = (A.dot( x_kmasuno ) - b)
            error.append( la.norm( difference, ord=np.inf))

        x_k = x_kmasuno
        tries -= 1

    if plot == True:
        domain = xrange(len(error))
        plt.semilogy(domain, error)
        plt.title("convergence of the Jacobi Method")
        plt.xlabel("Number of iterations", fontsize=14)
        plt.ylabel("Absolute value of approximation", fontsize=14)
        plt.show()

    return x_k

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
    A_a = np.copy(A)
    D = np.diag(A_a)
    d_inv= []
    for i in xrange(len(b)):
        d_inv.append(1./D[i])

    x_k = np.zeros(len(b))                      #changed to unsupported
    x_kmasuno = np.copy(x_k)
    this = False
    tries = maxiters
    error = []

    while this is False and tries > 0:

        for i in xrange(len(x_k)):
            x_kmasuno[i] = x_k[i] + d_inv[i]*(b[i] - A[i].T.dot(x_k))

        if ((la.norm((x_k - x_kmasuno), ord=np.inf)) <= tol):
            this = True
            difference = (A.dot( x_kmasuno ) - b)
            error.append( la.norm( difference, ord=np.inf))

        else:
            difference = (A.dot( x_kmasuno ) - b)
            error.append(la.norm( difference, ord=np.inf))

        x_k = np.copy(x_kmasuno)
        tries -= 1

    if plot == True:
        domain = xrange(len(error))
        plt.semilogy(domain, error)
        plt.title("convergence of the Gauss-Seidel")
        plt.xlabel("Number of iterations", fontsize=14)
        plt.ylabel("Absolute value of approximation", fontsize=14)
        plt.show()

    b = np.zeros_like((x_k))
    roar = np.column_stack((b,x_k))
    return roar[:,1:]


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
        A = diag_dom(i)
        b = np.random.rand(i,1)

        c = time.time()
        queso = la.solve(A,b)
        la_solve.append(time.time()-c)


        c = time.time()
        queso = gauss_seidel(A,b)
        gaus.append(time.time()-c)

    plt.loglog(domain, la_solve, 'g.-', linewidth=2, markersize=15, label = "lin solve")
    plt.loglog(domain, gaus, 'b.-', linewidth=2, markersize=15, label = "Gauss-Seidel")
    plt.title("Comparison", fontsize=18)
    plt.legend(loc="upper left")

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

    A_a = np.copy(A.toarray())
    D = np.diag(A_a)

    d_inv= []
    for i in xrange(len(b)):
        d_inv.append(1./D[i])

    x_k = np.zeros(len(b))                  #cambio de direccion
    x_kmasuno = np.copy(x_k)
    this = False
    tries = maxiters
    error = []

    while this is False and tries > 0:
        for i in xrange(len(x_k)):
            rowstart = A.indptr[i]
            rowend = A.indptr[i+1]
            Aix = np.dot(A.data[rowstart:rowend], x_k[A.indices[rowstart:rowend]])
            x_kmasuno[i] = x_k[i] + d_inv[i]*(b[i] - Aix)

        if ((la.norm((x_k - x_kmasuno), ord=np.inf)) <= tol):
            this = True
            difference = (A.dot( x_kmasuno ) - b)
            error.append( la.norm( difference, ord=np.inf))

        else:
            difference = (A.dot( x_kmasuno ) - b)
            error.append(la.norm( difference, ord=np.inf))

        x_k = np.copy(x_kmasuno)
        tries -= 1

    b = np.zeros_like((x_k))
    roar = np.column_stack((b,x_k))
    return roar[:,1:]


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

    x_kmasuno = np.zeros_like(b)
    x_k = np.copy(x_kmasuno)
    for j in xrange(maxiters):

        for i in xrange(len(x_k)):

            rowstart = A.indptr[i]
            rowend = A.indptr[i+1]
            Aix = np.dot(A.data[rowstart:rowend], x_kmasuno[A.indices[rowstart:rowend]])

            x_kmasuno[i] = x_k[i] + float(omega)*(1./A[i,i])*(b[i] - Aix)

        if la.norm((x_k - x_kmasuno), ord=np.inf) < tol:
            break

        x_k = np.copy(x_kmasuno)

    return x_kmasuno


# Problem 7
def finite_difference(n):
    """Return the A and b described in the finite difference problem that
    solves Laplace's equation.
    """
    offsets = [-1,0,1]
    B = sparse.diags([1,-4,1], offsets, shape=(n,n)).toarray()
    I = sparse.diags([1,1], [-n,n], shape=(n**2,n**2))
    A = []
    for i in xrange(n):
        A.append(B)
    uno = sparse.block_diag(A)

    b = np.zeros(n**2)
    b = b.astype(np.float64)
    b[::n] = -100
    b[n-1::n] = -100

    return (uno+I).tocsr(), b

# Problem 8
def compare_omega():
    """Time sparse_sor() with omega = 1, 1.05, 1.1, ..., 1.9, 1.95, tol=1e-2,
    and maxiters = 1000 using the A and b generated by finite_difference()
    with n = 20. Plot the times as a function of omega.
    """
    sparcify = []
    A,b = finite_difference(20)
    domain = np.arange(1,2,.05)
    for omega in domain:
        c = time.time()
        sparse_sor(A, b, omega, tol=1e-2, maxiters = 1000)
        D = time.time()
        sparcify.append(D-c)
        print omega

    plt.plot(domain, sparcify, 'g-', linewidth=2, markersize=15, label = "Timed SOR")
    plt.title("Timing", fontsize=18)
    plt.legend(loc="upper right")

    plt.show()


# Problem 9
def hot_plate(n):
    """Use finite_difference() to generate the system Au = b, then solve the
    system using SciPy's sparse system solver, scipy.sparse.linalg.spsolve().
    Visualize the solution using a heatmap using np.meshgrid() and
    plt.pcolormesh() ("seismic" is a good color map in this case).
    """
    #x = np.linspace(0, n, n)
    #y = x.copy()
    #X, Y = np.meshgrid(x, y)
    A,b = finite_difference(n)
    U = np.reshape(spla.spsolve(A,b),(n,n))
    plt.pcolormesh(U,cmap="seismic")
    plt.colorbar()

    plt.show()

if __name__ == "__main__":
    A = np.random.random((12,12))
    b = np.random.random(12)
    ros = gauss_seidel(A,b)
    print np.shape(ros)

    A_ = sparse.csr_matrix(diag_dom(5000))
    b_ = np.random.random(5000)
    roc = sparse_gauss_seidel(A_,b_)
    print np.shape(roc)


"""
    n=5
    #finite_difference(5)
    A = sparse.csr_matrix(diag_dom(5000))
    b = np.random.random(5000)
    a,b =  finite_difference(5)
    #print a.todense()
    #print b
    #print sparse_gauss_seidel(A,b)
    print compare_omega()

    # print prob4()
    #A = sparse.csr_matrix(diag_dom(20))
    #b = np.random.rand(20)
    #prob4()
    #a,b = finite_difference(4)
    #print a.toarray(),b
"""
