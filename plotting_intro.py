# matplotlib_intro.py
"""Introductory Labs: Intro to Matplotlib.
<Juan S Rodriguez>
<Math 345>
<September>
"""

import numpy as np
from matplotlib import pyplot as plt


def var_of_means(n):
    """Construct a random matrix A with values drawn from the standard normal
    distribution. Calculate the mean value of each row, then calculate the
    variance of these means. Return the variance.

    Inputs:
        n (int): The number of rows and columns in the matrix A.

    Returns:
        (float) The variance of the means of each row.
    """
    randmatrix = np.random.randn(n,n)
    mean = np.mean(randmatrix, axis=1)
    vari = np.var(mean)
    return vari


def prob1():
    """Create an array of the results of var_of_means() with inputs
    n = 100, 200, ..., 1000. Plot and show the resulting array.
    """
    var_array = []
    for i in xrange(100,1001,100):
        var_array.append(var_of_means(i))
    plt.plot(var_array)
    plt.show()
    return var_array


def prob2():
    """Plot the functions sin(x), cos(x), and arctan(x) on the domain
    [-2pi, 2pi]. Make sure the domain is refined enough to produce a figure
    with good resolution.
    """
    x = np.linspace(-2*(np.pi), 2*(np.pi), 100)
    y = np.cos(x)
    z = np.sin(x)
    w = np.arctan(x)
    plt.plot(x, y)
    plt.plot(x, z)
    plt.plot(x, w)
    plt.show()


def prob3():
    """Plot the curve f(x) = 1/(x-1) on the domain [-2,6].
        1. Split the domain so that the curve looks discontinuous.
        2. Plot both curves with a thick, dashed magenta line.
        3. Change the range of the y-axis to [-6,6].
    """
    x = np.linspace(-2, .9999, 50)
    y = 1/(x-1)
    plt.plot(x, y, 'm--', linewidth=5)
    x = np.linspace(1.0001, 6, 50)
    y = 1/(x-1)
    plt.plot(x, y, 'm--', linewidth=5)
    plt.ylim(-6,6)
    plt.show()


def prob4():
    """Plot the functions sin(x), sin(2x), 2sin(x), and 2sin(2x) on the
    domain [0, 2pi].
        1. Arrange the plots in a square grid of four subplots.
        2. Set the limits of each subplot to [0, 2pi]x[-2, 2].
        3. Give each subplot an appropriate title.
        4. Give the overall figure a title.
        5. Use the following line colors and styles.
              sin(x): green solid line.
             sin(2x): red dashed line.
             2sin(x): blue dashed line.
            2sin(2x): magenta dotted line.
    """
    plt.suptitle("Artistic Sinuses", fontsize=20)

    x = np.linspace(0, 2*(np.pi), 100)
    plt.subplot(221)
    plt.plot(x, np.sin(x), 'g-', lw=2)
    plt.axis([0, 2*(np.pi), -2, 2])
    plt.title("sin(x)", fontsize=12)

    plt.subplot(222)
    plt.plot(x, np.sin(2*x), 'r--', lw=2)
    plt.axis([0, 2*(np.pi), -2, 2])
    plt.title("sin(2x)", fontsize=12)

    plt.subplot(223)
    plt.plot(x, 2*np.sin(x), 'b--', lw=2)
    plt.axis([0, 2*(np.pi), -2, 2])
    plt.title("2sin(x)", fontsize=12)

    plt.subplot(224)
    plt.plot(x, 2*np.sin(2*x), 'm:', lw=2)
    plt.axis([0, 2*(np.pi), -2, 2])
    plt.title("2sin(2x)", fontsize=12)

    plt.show()


def prob5():
    """Visualize the data in FARS.npy. Use np.load() to load the data, then
    create a single figure with two subplots:
        1. A scatter plot of longitudes against latitudes. Because of the
            large number of data points, use black pixel markers (use "k,"
            as the third argument to plt.plot()). Label both axes.
        2. A histogram of the hours of the day, with one bin per hour.
            Label and set the limits of the x-axis.
    """
    fars = np.load("FARS.npy")
    plt.subplot(121)
    plt.plot(fars[:,1],fars[:,2],"k,")
    plt.axis("equal")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")

    plt.subplot(122)
    plt.hist(fars[:,0],  bins=np.arange(0, 23))
    plt.xlabel("Hours of the day")

    plt.show()


def prob6():
    """Plot the function f(x,y) = sin(x)sin(y)/xy on the domain
    [-2pi, 2pi]x[-2pi, 2pi].
        1. Create 2 subplots: one with a heat map of f, and one with a contour
            map of f. Choose an appropriate number of level curves, or specify
            the curves yourself.
        2. Set the limits of each subplot to [-2pi, 2pi]x[-2pi, 2pi].
        3. Choose a non-default color scheme.
        4. Add a colorbar to each subplot.
    """
    x = np.linspace(-2*np.pi, 2*np.pi, 200)
    y = x.copy()
    X, Y = np.meshgrid(x, y)

    Z = (np.sin(X) * np.sin(Y))/(X*Y)

    plt.subplot(121)
    plt.pcolormesh(X, Y, Z, cmap="viridis")
    plt.colorbar()
    plt.xlim(-2*np.pi, 2*np.pi)
    plt.ylim(-2*np.pi, 2*np.pi)

    plt.subplot(122)
    plt.contour(X, Y, Z, 30, cmap= "magma")
    plt.colorbar()
    plt.xlim(-2*np.pi, 2*np.pi)
    plt.ylim(-2*np.pi, 2*np.pi)

    plt.show()
