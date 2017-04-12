# qr_decomposition.py
"""Volume 1A: QR 1 (Decomposition).
<Juan>
<M345GL>
<October the last>
"""

import numpy as np
import scipy
from scipy import linalg as la


# Problem 1
def qr_gram_schmidt(A):
    """Compute the reduced QR decomposition of A via Modified Gram-Schmidt.

    Inputs:
        A ((m,n) ndarray): A matrix of rank n.

    Returns:
        Q ((m,n) ndarray): An orthonormal matrix.
        R ((n,n) ndarray): An upper triangular matrix.
    """
    m,n = np.shape(A)
    Q = np.copy(A)
    R = np.zeros((n,n))
    for i in xrange(n):
        R[i,i] = scipy.linalg.norm(Q[:,i])
        Q[:,i] = Q[:,i]/R[i,i]
        for j in xrange(i+1,n):
            R[i,j] = Q[:,j].T.dot(Q[:,i])
            Q[:,j] = Q[:,j] - R[i,j]*(Q[:,i])
    return Q,R

# Problem 2
def abs_det(A):
    """Use the QR decomposition to efficiently compute the absolute value of
    the determinant of A.

    Inputs:
        A ((n,n) ndarray): A square matrix.

    Returns:
        (float) the absolute value of the detetminant of A.
    """
    return np.absolute(np.prod(np.diag(qr_gram_schmidt(A)[1])))

# Problem 3
def solve(A, b):
    """Use the QR decomposition to efficiently solve the system Ax = b.

    Inputs:
        A ((n,n) ndarray): An invertible matrix.
        b ((n, ) ndarray): A vector of length n.

    Returns:
        x ((n, ) ndarray): The solution to the system Ax = b.
    """
    Q,R = la.qr(A, mode="economic")
    Y = Q.T.dot(b)

    m = len(Y)
    X = [0]*(m)
    for k in xrange(m-1,-1,-1):             #k represents the row
        sumb = 0
        for l in xrange(m-1, k, -1):         #l represents the column
            sumb += R[k,l]*X[l]
        X[k]=(Y[k]-sumb)*(1/R[k][k])
    return X

# Problem 4
def qr_householder(A):
    """Compute the full QR decomposition of A via Householder reflections.

    Inputs:
        A ((m,n) ndarray): A matrix of rank n.

    Returns:
        Q ((m,m) ndarray): An orthonormal matrix.
        R ((m,n) ndarray): An upper triangular matrix.
    """
    sign = lambda x: 1 if x >= 0 else -1

    m,n = np.shape(A)
    R = np.copy(A)
    Q = np.identity(m)

    for k in xrange(n):

        u = np.copy(R[k:,k])
        u[0] = u[0] + sign(u[0])*scipy.linalg.norm(u)
        u = u/scipy.linalg.norm(u)
        R[k:,k:] = R[k:,k:] - 2*np.outer(u,(u.T.dot(R[k:,k:])))
        Q[k:,:] = Q[k:,:] - 2*np.outer(u,(u.T.dot(Q[k:,:])))

    return Q.T, R

# Problem 5
def hessenberg(A):
    """Compute the Hessenberg form H of A, along with the orthonormal matrix Q
    such that A = QHQ^T.

    Inputs:
        A ((n,n) ndarray): An invertible matrix.

    Returns:
        H ((n,n) ndarray): The upper Hessenberg form of A.
        Q ((n,n) ndarray): An orthonormal matrix.
    """
    m,n = np.shape(A)
    H = np.copy(A)
    Q = np.identity(n)
    for k in xrange(n-2):
        u = np.copy(H[(k+1):, k])
        u[0] = u[0] + scipy.linalg.norm(u)
        u = u/scipy.linalg.norm(u)

        H[(k+1):, k:] = H[(k+1):, k:] - np.outer(2*u,(np.dot(np.transpose(u),(H[(k+1):, k:]))))
        H[:,k+1:] = H[:,k+1:] - 2*np.outer((H[:,(k+1):].dot(u)), (u.T))
        Q[k+1:,:] = Q[k+1:,:] - 2*np.outer(u,(u.T.dot(Q[(k+1):,:])))

    return H, np.transpose(Q)

    """

if __name__ == "__main__":
    #prob4
    A = np.random.random((5, 3))
    Q,R = la.qr(A)
    q,r = qr_householder(A)
    print A.shape, q.shape, r.shape
    print np.allclose(q.dot(r), A)

    #prob5
    A = np.random.random((8,8))
    H, Q = la.hessenberg(A, calc_q=True)
    h, q = hessenberg(A)

    print np.allclose(np.triu(h, -1), h)
    print np.allclose(np.dot(np.dot(q, h), q.T), A)


    #Problem 1 tester
    B = np.copy(A)
    Q,R = la.qr(A, mode="economic")
    S,T = qr_gram_schmidt(B)
    print Q
    print S
    print R

    print T
    print np.allclose(np.dot(S, T), B)

    #problem 3 tester
    A_ = np.copy(A)
    b_ = np.copy(b)
    print solve(A,b)
    print la.solve(A_,b_)
    """
