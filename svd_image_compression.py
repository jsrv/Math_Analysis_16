# svd_image_compression.py
"""Volume 1A: SVD and Image Compression.
<Name>
<Class>
<Date>
"""

from scipy import linalg as la
from numpy import linalg as lak
import numpy as np
from matplotlib import pyplot as plt

# Problem 1
def truncated_svd(A,k=None):
    """Computes the truncated SVD of A. If r is None or equals the number
        of nonzero singular values, it is the compact SVD.
    Parameters:
        A: the matrix
        k: the number of singular values to use
    Returns:
        U - the matrix U in the SVD
        s - the diagonals of Sigma in the SVD
        Vh - the matrix V^H in the SVD
    """
    U, s, V_h = la.svd(A)
    #print U
    #print s
    #print V_h
    rk = lak.matrix_rank(A)
    m,n = A.shape
    A_h = np.transpose(np.conj(A))
    AHA = A_h.dot(A)
    eigval, eigvec = la.eig(AHA)

    idx = eigval.argsort()[::-1]
    eigval = eigval[idx]
    eigvec = eigvec[:,idx]

    sing_val = np.sqrt(eigval)
    sing_val = np.trim_zeros(sing_val)
    if k != None:
        sing_val = sing_val[:k]
    V = eigvec[:,:len(sing_val)]
    U = np.dot(A,V)
    m,n = np.shape(U)

    for i in xrange(n):
        U[:,i] = U[:,i]*(1./sing_val[i])

    return U, sing_val, V


# Problem 2
def visualize_svd():
    """Plot each transformation associated with the SVD of A."""
    A = np.array([[3,1],[1,3]])
    U, sig, V_h = la.svd(A)
    s =np.diag(sig)
    plt.suptitle("Transformations by SVD", fontsize=20)

    e_1 = np.array([1,0])
    e_2 = np.array([0,1])

    x = np.linspace(0, 2*(np.pi), 100)
    exis = np.cos(x)
    yes = np.sin(x)
    this = np.vstack((exis,yes))
    otra = V_h.dot(this)
    e_1 = np.array([[0,0],[0,1]])
    e_2 = np.array([[0,1],[0,0]])
    vectors1 = V_h.dot(e_1)
    vectors2 = V_h.dot(e_2)

    plt.subplot(221)
    plt.plot(exis, yes, 'b-', lw=1)
    plt.plot(e_1[0,:], e_1[1,:], 'g-', lw=1)
    plt.plot(e_2[0,:], e_2[1,:], 'g-', lw=1)
    plt.axis("equal")
    plt.title("S", fontsize=12)

    plt.subplot(222)
    plt.plot(otra[0,:], otra[1,:], 'b-', lw=1)
    plt.plot(vectors1[0,:], vectors1[1,:], 'g-', lw=1)
    plt.plot(vectors2[0,:], vectors2[1,:], 'g-', lw=1)
    plt.axis("equal")
    plt.title("VhS", fontsize=12)

    plt.subplot(223)
    plt.plot(s.dot(otra)[0,:], s.dot(otra)[1,:], 'b-', lw=1)
    plt.plot(s.dot(V_h).dot(e_1)[0,:], s.dot(V_h).dot(e_1)[1,:], 'g-', lw=1)
    plt.plot(s.dot(V_h).dot(e_2)[0,:], s.dot(V_h).dot(e_2)[1,:], 'g-', lw=1)
    plt.axis("equal")
    plt.title("ZVhS", fontsize=12)

    plt.subplot(224)
    plt.plot(U.dot(s).dot(otra)[0,:], U.dot(s).dot(otra)[1,:], 'b-', lw=1)
    plt.plot(U.dot(s).dot(V_h).dot(e_1)[0,:], U.dot(s).dot(V_h).dot(e_1)[1,:], 'g-', lw=1)
    plt.plot(U.dot(s).dot(V_h).dot(e_2)[0,:], U.dot(s).dot(V_h).dot(e_2)[1,:], 'g-', lw=1)
    plt.axis("equal")
    plt.title("UZVh", fontsize=12)

    plt.axis("equal")
    plt.show()

# Problem 3
def svd_approx(A, k):
    """Returns best rank k approximation to A with respect to the induced 2-norm.

    Inputs:
    A - np.ndarray of size mxn
    k - rank

    Return:
    Ahat - the best rank k approximation
    """
    U, sig, V_h = la.svd(A, full_matrices=False)
    S = np.diag(sig[:k])
    Ahat = U[:,:k].dot(S).dot(V_h[:k,:])
    return Ahat


# Problem 4
def lowest_rank_approx(A,e):
    """Returns the lowest rank approximation of A with error less than e
    with respect to the induced 2-norm.

    Inputs:
    A - np.ndarray of size mxn
    e - error

    Return:
    Ahat - the lowest rank approximation of A with error less than e.
    """
    U, sig, V_h = la.svd(A, full_matrices=False)
    print sig
    rk = lak.matrix_rank(A)
    Ahat = A
    dif = 0

    while dif < e:
        dif = la.norm(A-Ahat)
        if dif >= e:
            break
        car = np.copy(Ahat)
        Ahat = U[:,:rk].dot(np.diag(sig[:rk])).dot(V_h[:rk,:])
        rk -= 1
        if rk <= 0:
            break
    return car


# Problem 5
def compress_image(filename,k):
    """Plot the original image found at 'filename' and the rank k approximation
    of the image found at 'filename.'

    filename - jpg image file path
    k - rank
    """
    my_file = plt.imread(filename)[:,:,:].astype(float)

    R = my_file[:,:,0]
    G = my_file[:,:,1]
    B = my_file[:,:,2]
    Krazo = np.zeros_like(my_file)

    Krazo[:,:,0] = svd_approx(R,k)
    Krazo[:,:,1] = svd_approx(G,k)
    Krazo[:,:,2] = svd_approx(B,k)

    my_file = np.round(my_file)/255.
    K = np.round(Krazo)/255.
    K[K<0] = 0
    K[K>1] = 1
    plt.subplot(121)
    plt.imshow(my_file)
    plt.subplot(122)
    plt.imshow(K)
    plt.show()

if __name__ == "__main__":
    a = np.array([[ -8. , 7. ,  0. , -6. ,  2. , -6. ,  5. ,  8.,   1.],
 [ -5.,  -8.,   4.,  -8.,  -7.,  -2.,   5.,   8. ,  4.],
 [ -3. ,  0.,  -2.,   9.,   6.,   9.,   3.,  -7.,  -5.],
 [ -9.  , 6.,   1., -10.,  -6.,  -7.,   8.,  -7.,   2.],
 [ -3.  ,-1.,   8.,  -2.,   3.,  -8.,  -3., -10.,   8.],
 [ -3.,   7.,   0.,  -7.,  -4.,  -1.,   2.,  -8.,   4.]])
    b = np.array([[1,0,0,0,2],[0,0,3,0,0],[0,0,0,0,0],[0,2,0,0,0]])
    U, s, V = truncated_svd(b)
    S = np.diag(s)
    Aprox = U.dot(S.dot(V.T))
    
    print np.allclose(Aprox,b)

"""
    #print visualize_svd()
    compress_image('hubble_image.jpg',20)

    a = np.array([[1,0,0,0,2],[0,0,3,0,0],[0,0,0,0,0],[0,2,0,0,0]])
    U, s, V = truncated_svd(a)
    S = np.diag(s)
    Aprox = U.dot(S.dot(V.T))
    d = np.copy(a)
    print Aprox
    print d
    print np.allclose(Aprox, d)

    d = np.copy(a)
    U,s,V = truncated_svd(a)
    print
    print U
    print s
    print V
    print

    S = np.diag(s)
    Aprox = U.dot(S.dot(V.T))
    print Aprox
    print d
    print np.allclose(Aprox, d)
"""
