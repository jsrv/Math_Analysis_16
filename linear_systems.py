# linear_systems.py
"""Volume 1A: Linear Systems.
<Juan>
<M345G>
<Oct/2016>
"""
import numpy as np
import random
import time
import scipy
from numpy import array
from matplotlib import pyplot as plt
from scipy import sparse
from scipy import linalg
from scipy.sparse import linalg as spla

# Problem 1
def ref(A):
    """Reduce the square matrix A to REF. You may assume that A is invertible
    and that a 0 will never appear on the main diagonal. Avoid operating on
    entries that you know will be 0 before and after a row operation.
    """
    m, n = A.shape
    for i in (xrange(m)):
        for j in (xrange(i+1,m)):
            #j as row and i as column
            A[j,i:] -= (A[j,i] / A[i,i]) * A[i,i:]
            #A[2,1:] -= (A[2,1] / A[1,1]) * A[1,1:]
    return A


# Problem 2
def lu(A):
    """Compute the LU decomposition of the square matrix A. You may assume the
    decomposition exists and requires no row swaps.

    Returns:
        L ((n,n) ndarray): The lower-triangular part of the decomposition.
        U ((n,n) ndarray): The upper-triangular part of the decomposition.
    """
    m, n = A.shape
    U = np.copy(A)
    L = np.identity(m)
    for j in (xrange(n)):
        for i in (xrange(j+1,m)):
            L[i,j] = U[i][j]/U[j][j]
            #U[j,i:] -= (A[j,i] / A[i,i]) * A[i,i:]
            U[i,j:] -= L[i,j]*U[j,j:]
    return L,U


# Problem 3
def solve(A, b):
    """Use the LU decomposition and back substitution to solve the linear
    system Ax = b. You may assume that A is invertible (hence square).
    """
    Y = []
    m, n = A.shape
    L, U = lu(A)
    for i in xrange(m):             #i represents the row
        suma = 0

        for j in xrange(i):         #j represents the column
            suma += L[i,j]*Y[j]
        Y.append(b[i]-suma)


    m = len(Y)
    X = [0]*(m)
    for k in xrange(m-1,-1,-1):             #k represents the row
        sumb = 0
        for l in xrange(m-1, k, -1):         #l represents the column
            print k,l
            sumb += U[k,l]*X[l]
        X[k]=(Y[k]-sumb)*(1/U[k][k])
    return X


# Problem 4
def prob4():
    """Time different scipy.linalg functions for solving square linear systems.
    Plot the system size versus the execution times. Use log scales if needed.
    """
    domain = 2**np.arange(1,12)
    lainv = []
    lasol = []
    lalu1 = []
    lalu2 = []
    for i in domain:
        A = 5*np.random.rand(i,i)
        b = 5*np.random.rand(i)

        Tsolve = time.time()
        inv = scipy.linalg.inv(A)
        sol = inv.dot(b)
        lainv.append(time.time()-Tsolve)

        Tsolve = time.time()
        sol = scipy.linalg.solve(A,b)
        lasol.append(time.time()-Tsolve)

        Tsolve = time.time()
        lu, piv = linalg.lu_factor(A)
        linalg.lu_solve((lu,piv), b)
        lalu1.append(time.time()-Tsolve)


        lu, piv = linalg.lu_factor(A)
        Tsolve = time.time()
        linalg.lu_solve((lu,piv), b)
        lalu2.append(time.time()-Tsolve)


    plt.subplot(221)
    plt.plot(domain, lainv)
    plt.title("Inverse", fontsize=18)

    plt.subplot(222)
    plt.plot(domain, lasol)
    plt.title("Solve", fontsize=18)

    plt.subplot(223)
    plt.plot(domain, lalu1)
    plt.title("LU Factor and Solve", fontsize=18)

    plt.subplot(224)
    plt.plot(domain, lalu2)
    plt.title("LU Solve", fontsize=18)
    plt.show()


# Problem 5
def prob5(n):
    """Return a sparse n x n tridiagonal matrix with 2's along the main
    diagonal and -1's along the first sub- and super-diagonals.
    """
    return scipy.sparse.diags([-1,2,-1], [-1,0,1], shape=(n,n)).toarray()


# Problem 6
def prob6():
    """Time regular and sparse linear system solvers. Plot the system size
    versus the execution times. As always, use log scales where appropriate.
    """
    domain = 2**np.arange(1,13)
    spaser = []
    denser = []
    for i in domain:
        b = 5*np.random.rand(i)
        A = prob5(i)
        B = np.copy(A)
        Acsr = scipy.sparse.csr_matrix(A)

        Tsolve = time.time()
        R = scipy.sparse.linalg.spsolve(Acsr,b)
        spaser.append(time.time()-Tsolve)

        Tsolve = time.time()
        S = scipy.linalg.solve(B,b)
        denser.append(time.time()-Tsolve)

    plt.subplot(121)
    plt.plot(domain, spaser)
    plt.title("Sparsed Solve", fontsize=18)

    plt.subplot(122)
    plt.plot(domain, denser)
    plt.title("Dense Solve", fontsize=18)

    plt.show()
