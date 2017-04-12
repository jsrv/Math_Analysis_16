# linear_transformations.py
"""Volume 1A: Linear Transformations.
<Juan Sebastian Rodriguez>
<M345L>
<September 20, 2016>
"""

import numpy as np
from random import random
from matplotlib import pyplot as plt
import time


# Problem 1
def notprob(orig_imag,change_image):
    plt.subplot(121)
    plt.title("original image")
    plt.plot(orig_imag[0], orig_imag[1], 'k,')
    plt.axis([-1,1,-1,1])
    plt.gca().set_aspect("equal")

    plt.subplot(122)
    plt.title("Changed image")
    plt.plot(change_image[0], change_image[1], 'k,')
    plt.axis([-1,1,-1,1])
    plt.gca().set_aspect("equal")

    plt.show()



def stretch(A, a, b):
    """Scale the points in 'A' by 'a' in the x direction and 'b' in the
    y direction.

    Inputs:
        A ((2,n) ndarray): Array containing points in R2 stored as columns.
        a (float): scaling factor in the x direction.
        b (float): scaling factor in the y direction.
    """
    change = np.array([[a,0],[0,b]])
    trans_array = np.dot(change,A)
    return trans_array

def shear(A, a, b):
    """Slant the points in 'A' by 'a' in the x direction and 'b' in the
    y direction.

    Inputs:
        A ((2,n) ndarray): Array containing points in R2 stored as columns.
        a (float): scaling factor in the x direction.
        b (float): scaling factor in the y direction.
    """
    change = np.array([[1,a],[b,1]])
    trans_array = np.dot(change,A)
    return trans_array

def reflect(A, a, b):
    """Reflect the points in 'A' about the line that passes through the origin
    and the point (a,b).

    Inputs:
        A ((2,n) ndarray): Array containing points in R2 stored as columns.
        a (float): x-coordinate of a point on the reflecting line.
        b (float): y-coordinate of the same point on the reflecting line.
    """
    change = np.array([[a**2-b**2,2*a*b],[2*a*b,b**2-a**2]])
    trans_array = np.dot(change,A)
    trans_arraydiv = trans_array/float(a**2+b**2)
    return trans_arraydiv

def rotate(A, theta):
    """Rotate the points in 'A' about the origin by 'theta' radians.

    Inputs:
        A ((2,n) ndarray): Array containing points in R2 stored as columns.
        theta (float): The rotation angle in radians.
    """
    change = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    trans_array = np.dot(change,A)
    return trans_array


# Problem 2
def solar_system(T, omega_e, omega_m):
    """Plot the trajectories of the earth and moon over the time interval [0,T]
    assuming the initial position of the earth is (10,0) and the initial
    position of the moon is (11,0).

    Parameters:
        T (int): The final time.
        omega_e (float): The earth's angular velocity.
        omega_m (float): The moon's angular velocity.
    """
    Pe0 = np.array([10,0])
    Pm0 = np.array([11,0])
    t = np.linspace(0,T, 500)
    Petar = []
    Pmtar = []

    for i in t:

        Pet = rotate(Pe0,i*omega_e)
        Petar.append(Pet)

        Pmrot = rotate((Pm0-Pe0), i*omega_m)
        Pmtar.append(Pmrot + Petar[-1])

    earthx, earthy = zip(*Petar)
    moonx, moony = zip(*Pmtar)
    plt.plot(earthx, earthy, label='Earth bro')
    plt.plot(moonx, moony, label='Moon bro')
    plt.legend(loc="upper left")
    plt.gca().set_aspect("equal")
    plt.show()



def random_vector(n):
    """Generate a random vector of length n as a list."""
    return [random() for i in xrange(n)]

def random_matrix(n):
    """Generate a random nxn matrix as a list of lists."""
    return [[random() for j in xrange(n)] for i in xrange(n)]

def matrix_vector_product(A, x):
    """Compute the matrix-vector product Ax as a list."""
    m, n = len(A), len(x)
    return [sum([A[i][k] * x[k] for k in range(n)]) for i in range(m)]

def matrix_matrix_product(A, B):
    """Compute the matrix-matrix product AB as a list of lists."""
    m, n, p = len(A), len(B), len(B[0])
    return [[sum([A[i][k] * B[k][j] for k in range(n)])
                                    for j in range(p) ]
                                    for i in range(m) ]


# Problem 3
def prob3():
    """Use time.time(), timeit.timeit(), or %timeit to time
    matrix_vector_product() and matrix-matrix-mult() with increasingly large
    inputs. Generate the inputs A, x, and B with random_matrix() and
    random_vector() (so each input will be nxn or nx1).
    Only time the multiplication functions, not the generating functions.

    Report your findings in a single figure with two subplots: one with matrix-
    vector times, and one with matrix-matrix times. Choose a domain for n so
    that your figure accurately describes the growth, but avoid values of n
    that lead to execution times of more than 1 minute.
    """
    domain = 2**np.arange(1,9)
    times1 = []
    times2 = []
    times3 = []
    for n in domain:
        randmat = random_matrix(n)
        randmatB = random_matrix(n)
        rand_vec = random_vector(n)

        start1 = time.time()
        matrix_vector_product(randmat,rand_vec)
        times1.append(time.time() - start1)

        start2 = time.time()
        matrix_matrix_product(randmat,randmatB)
        times2.append(time.time() - start2)

    plt.subplot(121)
    plt.title("Matrix Vector multiplication")
    plt.plot(domain, times1, 'g.-', linewidth=2, markersize=15)

    plt.subplot(122)
    plt.title("Matrix Matrix multiplication")
    plt.plot(domain, times2, 'g.-', linewidth=2, markersize=15)
    plt.xlabel("n", fontsize=14)
    plt.ylabel("Seconds", fontsize=14)
    plt.show()


# Problem 4
def prob4(N=8):
    """Time matrix_vector_product(), matrix_matrix_product(), and np.dot().

    Report your findings in a single figure with two subplots: one with all
    four sets of execution times on a regular linear scale, and one with all
    four sets of exections times on a log-log scale.
    """
    domain = 2**np.arange(1,N+1)
    times1 = []
    times2 = []
    times3 = []
    times4 = []
    for n in domain:
        randmat = random_matrix(n)
        randmatB = random_matrix(n)
        rand_vec = random_vector(n)

        start1 = time.time()
        matrix_vector_product(randmat,rand_vec)
        times1.append(time.time() - start1)

        start2 = time.time()
        matrix_matrix_product(randmat,randmatB)
        times2.append(time.time() - start2)

        start3 = time.time()
        np.dot(randmat, rand_vec)
        times3.append(time.time() - start3)

        start4 = time.time()
        np.dot(randmat, randmatB)
        times4.append(time.time() - start4)

    plt.subplot(121)
    plt.title("Normal Scale")
    plt.plot(domain, times1, 'b.-', linewidth=2, markersize=15)
    plt.plot(domain, times2, 'g.-', linewidth=2, markersize=15)
    plt.plot(domain, times3, 'r.-', linewidth=2, markersize=15)
    plt.plot(domain, times4, 'y.-', linewidth=2, markersize=15)

    plt.subplot(122)
    plt.title("log-log scale")
    plt.loglog(domain, times1, 'b.-', linewidth=2, markersize=15)
    plt.loglog(domain, times2, 'g.-', linewidth=2, markersize=15)
    plt.loglog(domain, times3, 'r.-', linewidth=2, markersize=15)
    plt.loglog(domain, times4, 'y.-', linewidth=2, markersize=15)


    plt.xlabel("n", fontsize=14)
    plt.ylabel("Seconds", fontsize=14)
    plt.show()
