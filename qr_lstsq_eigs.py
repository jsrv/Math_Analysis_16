# qr_lstsq_eigs.py
"""Volume 1A: QR 2 (Least Squares and Computing Eigenvalues).
<Name>
<Class>
<Date>"""

import numpy as np
import cmath
from scipy import linalg as la
from matplotlib import pyplot as plt
from numpy import poly1d as pol
import random



# Problem 1
def least_squares(A, b):
    """Calculate the least squares solutions to Ax = b using QR decomposition.

    Inputs:
        A ((m,n) ndarray): A matrix of rank n <= m.
        b ((m, ) ndarray): A vector of length m.

    Returns:
        x ((n, ) ndarray): The solution to the normal equation.
    """
    Q, R = la.qr(A, mode = "economic")
    x_hat = la.solve_triangular(R,(np.dot(np.transpose(Q), b)))
    return x_hat



# Problem 2
def line_fit():
    """Load the data from housing.npy. Use least squares to calculate the line
    that best relates height to weight.

    Plot the original data points and the least squares line together.
    """
    data = np.load("housing.npy")
    year, index = np.load("housing.npy").T
    b = data[:,1:]
    D = np.copy(data)

    A = np.column_stack((D[:,0], np.ones_like(D[:,0])))
    least_squares_sol = least_squares(A,b)
    plt.plot(year, index, "o", markersize = 8)
    plt.plot(year, A.dot(least_squares_sol))
    plt.show()


# Problem 3
def polynomial_fit():
    """Load the data from housing.npy. Use least squares to calculate
    the polynomials of degree 3, 6, 9, and 12 that best fit the data.

    Plot the original data points and each least squares polynomial together
    in individual subplots.
    """
    data = np.load("housing.npy")
    year, index = np.load("housing.npy").T
    arco = []
    b = data[:,1:]
    D = np.copy(data)
    A_ = np.column_stack((D[:,0], np.ones_like(D[:,0])))
    this_ = D[:,0]
    for i in xrange(11):
        arco.append(A_)
        this_mas = np.multiply(this_, D[:,0])
        A_mas = np.insert(A_, 0, this_mas, axis=1)
        A_ = A_mas
        this_ = this_mas
    arco.append(A_)
    a = 0
    for i in ([2,5,8,11]):
        a += 1
        plt.subplot(2,2,a)
        x = la.lstsq(arco[i], b)[0]
        plt.plot(year, arco[i].dot(x))
        plt.plot(year, index, "o", markersize = 3)
    x = la.lstsq(arco[3], b)[0]
    plt.show()



def plot_ellipse(a, b, c, d, e):
    """Plot an ellipse of the form ax^2 + bx + cxy + dy + ey^2 = 1."""
    xk, yk = np.load("ellipse.npy").T


    theta = np.linspace(0, 2*np.pi, 200)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    A = a*(cos_t**2) + c*cos_t*sin_t + e*(sin_t**2)
    B = b*cos_t + d*sin_t
    r = (-B + np.sqrt(B**2 + 4*A))/(2*A)

    plt.plot(r*cos_t, r*sin_t, lw=2)
    plt.gca().set_aspect("equal", "datalim")
    plt.plot(xk, yk, "o")
    plt.show()



# Problem 4
def ellipse_fit():
    """Load the data from ellipse.npy. Use least squares to calculate the
    ellipse that best fits the data.

    Plot the original data points and the least squares ellipse together.
    """
    xk, yk = np.load("ellipse.npy").T
    A = np.column_stack((xk**2, xk, np.multiply(xk,yk), yk, yk**2))
    b = np.ones_like(xk)
    a, b, c, d, e = la.lstsq(A, b)[0]
    plot_ellipse(a, b, c, d, e)

# Problem 5
def power_method(A, N=20, tol=1e-12):
    """Compute the dominant eigenvalue of A and a corresponding eigenvector
    via the power method.

    Inputs:
        A ((n,n) ndarray): A square matrix.
        N (int): The maximum number of iterations.
        tol (float): The stopping tolerance.

    Returns:
        (foat): The dominant eigenvalue of A.
        ((n, ) ndarray): An eigenvector corresponding to the dominant
            eigenvalue of A.
    """
    eigs, vecs = la.eig(A)
    loc = np.argmax(eigs)
    lamb, x = eigs[loc], vecs[:,loc]
    #print lamb
    #print x

    m, n = np.shape(A)
    x_ = np.random.rand(n)
    x_ = x_/la.norm(x_)
    for k in xrange(N):
        x_1 = A.dot(x_)
        x_1 = x_1/la.norm(x_1)

        if la.norm((x_1 - x_)) <= tol:
            break

        x_ = x_1.copy()

    return x_.dot(A.dot(x_))/x_.dot(x_), x_1


# Problem 6
def qr_algorithm(A, N=50, tol=1e-12):
    """Compute the eigenvalues of A via the QR algorithm.

    Inputs:
        A ((n,n) ndarray): A square matrix.
        N (int): The number of iterations to run the QR algorithm.
        tol (float): The threshold value for determining if a diagonal block
            is 1x1 or 2x2.

    Returns:
        ((n, ) ndarray): The eigenvalues of A.
    """
    m,n = np.shape(A)
    S = la.hessenberg(A)
    for k in xrange(N):
        Q,R = la.qr(S, mode = "economic")
        S = R.dot(Q)
    n,n = np.shape(S)
    eigs = []
    i = 0
    while i < n:
        if ([i+1,i+1]==[n,n]) or (np.abs(S[i+1,i]) < tol) :
            eigs.append(S[i,i])
        else:
            mat = np.array([[S[i,i],S[i,i+1]],[S[i+1,i],S[i+1,i+1]]])
            eigens = la.eig(mat)[0]
            #print eigens
            #b = (-S[i,i])*S[i+1,i+1]
            #c = (S[i,i]*S[i+1,i+1] - S[i,i+1]*S[i+1,i])
            #r = (-b + cmath.sqrt(b**2 - 4*c)) / (2.)
            #s = (-b - cmath.sqrt(b**2 - 4*c)) / (2.)
            #print
            #print r,s
            #print eigens
            eigs.append(float(eigens[1]))
            eigs.append(float(eigens[0]))

            i += 1
        i += 1
    print len(eigs)
    return eigs

"""
if __name__ == "__main__":
    A = np.random.random((10,10))
    b = A.T + A
    eigens = la.eig(b)[0]
    print list(eigens)
    print qr_algorithm(b)
"""
