# numpy_intro.py
"""Introductory Labs: Intro to NumPy.
<Juan>
<345>
<Today>
"""
import numpy as np

def prob1():
    """Define the matrices A and B as arrays. Return the matrix product AB."""
    A = np.array([[3,-1,4],[1,5,-9]])
    B = np.array([[2,6,-5,3],[5,-8,9,7],[9,-3,-2,-3]])
    C = np.dot(A,B)
    return C


def prob2():
    """Define the matrix A as an array. Return the matrix -A^3 + 9A^2 - 15A."""
    A = np.array([[3,1,4],[1,5,9],[-5,3,1]])
    onedot = np.dot(A,A)
    result = -(np.dot(A,onedot)) + 9*(onedot) - 15*A
    return result


def prob3():
    """Define the matrices A and B as arrays. Calculate the matrix product ABA,
    change its data type to np.int64, and return it.
    """
    D = np.ones((7,7),dtype=np.int)
    A = np.triu(D)
    B = 5*D - 6*np.tril(D)
    print D
    print A
    print B
    product = A.dot(B)
    product2 = product.dot(A)
    pprod = product2.astype(np.int64)
    return pprod


def prob4(aray):
    """Make a copy of 'A' and set all negative entries of the copy to 0.
    Return the copy.

    Example:
        >>> A = np.array([-3,-1,3])
        >>> prob4(A)
        array([0, 0, 3])
    """
    B = np.array(aray)
#    B = A
    mask = B <= 0
    B[mask] = 0
    return B


def prob5():
    """Define the matrices A, B, and C as arrays. Return the block matrix
                                | 0 A^T I |
                                | A  0  0 |,
                                | B  0  C |
    where I is the identity matrix of appropriate size and each 0 is a matrix
    of all zeros, also of appropriate sizes.
    """
    zeros2x2= np.zeros(4)
    zeros2x2 = zeros2x2.reshape((2,2))
    zeros3x3= np.zeros(9)
    zeros3x3 = zeros3x3.reshape((3,3))
    zeros3x2= np.zeros(6)
    zeros3x2 = zeros3x2.reshape((3,2))
    ident = np.eye(3)
    A = np.array([[0,2,4],[1,3,5]])
    B = np.array([[3,0,0],[3,3,0],[3,3,3]])
    C = (-2)*np.eye(3)
    AT = A.T
    linea = np.hstack((zeros3x3, AT, ident))
    lineb = np.hstack((A, zeros2x2, zeros3x2.T))
    linec = np.hstack((B, zeros3x2, C))
    matrixs = np.vstack((linea, lineb, linec))
    return matrixs


def prob6(A):
    """Divide each row of 'A' by the row sum and return the resulting array.

    Example:
        >>> A = np.array([[1,1,0],[0,1,0],[1,1,1]])
        >>> prob6(A)
        array([[ 0.5       ,  0.5       ,  0.        ],
               [ 0.        ,  1.        ,  0.        ],
               [ 0.33333333,  0.33333333,  0.33333333]])
    """
    A = A.astype(float)
    result = np.divide(A,A.sum(axis=1)[:,None])
    return result


def prob7():
    """Given the array stored in grid.npy, return the greatest product of four
    adjacent numbers in the same direction (up, down, left, right, or
    diagonally) in the grid.
    """
    grid = np.load("grid.npy")
    a=(np.max(grid[:,:-3] * grid[:,1:-2] * grid[:,2:-1] * grid[:,3:]))
    b=(np.max(grid[:-3,:] * grid[1:-2,:] * grid[2:-1,:] * grid[3:,:]))
    c=(np.max(grid[:-3,:-3] * grid[1:-2,1:-2] * grid[2:-1,2:-1] * grid[3:,3:]))
    d=(np.max(grid[:-3,3:] * grid[1:-2,2:-1] * grid[2:-1,1:-2] * grid[3:,:-3]))
    e = max(a,b,c,d)
    print grid
    return e


if __name__ == "__main__":
    print prob3()
