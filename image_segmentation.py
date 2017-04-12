# image_segmentation.py
"""Volume 1A: Image Segmentation.
<Name>
<Class>
<Date>
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy import linalg as la
from scipy.sparse.linalg import eigsh
from scipy.sparse import lil_matrix
from collections import Counter
from scipy import sparse as sp
import copy
import heapq


# Problem 1: Implement this function.
def laplacian(A):
    '''
    Compute the Laplacian matrix of the adjacency matrix A.
    Inputs:
        A (array): adjacency matrix for undirected weighted graph,
             shape (n,n)
    Returns:
        L (array): Laplacian matrix of A

    '''
    m,n = A.shape
    A = -1*(A)
    for row in xrange(m):
        suma = sum(A[:,row])
        A[row,row] = -1 * suma
    return A

# Problem 2: Implement this function.
def n_components(A,tol=1e-8):
    '''
    Compute the number of connected components in a graph
    and its algebraic connectivity, given its adjacency matrix.
    Inputs:
        A -- adjacency matrix for undirected weighted graph,
             shape (n,n)
        tol -- tolerance value
    Returns:
        n_components -- the number of connected components
        lambda -- the algebraic connectivity
    '''
    lap = laplacian(A)
    real_eval = np.real(la.eigh(lap)[0])
    real_eval[real_eval<tol] = 0

    zeros = np.count_nonzero(real_eval == 0)
    connectivity = 0
    if zeros == 1:
        minima = heapq.nsmallest(2,real_eval)
        connectivity = max(minima)

    return zeros, connectivity
    '''
    lap = laplacian(A)
    eigval, eigvec = la.eigh(lap)
    real_eval = np.real(eigval)
    copy_reale = list(real_eval)
    num_zeros = 0

    for i in xrange(len(real_eval)):
        if real_eval[i] < tol:
            num_zeros += 1
            copy_reale.remove(real_eval[i])

    minima = min(copy_reale)
    c = Counter(real_eval)
    return num_zeros, minima
    '''
# Problem 3: Implement this function.
def adjacency(filename="dream.png", radius = 5.0, sigma_I = .02, sigma_d = 3.0):
    '''
    Compute the weighted adjacency matrix for
    the image given the radius. Do all computations with sparse matrices.
    Also, return an array giving the main diagonal of the degree matrix.

    Inputs:
        filename (string): filename of the image for which the adjacency matrix will be calculated
        radius (float): maximum distance where the weight isn't 0
        sigma_I (float): some constant to help define the weight
        sigma_d (float): some constant to help define the weight
    Returns:
        W (sparse array(csc)): the weighted adjacency matrix of img_brightness,
            in sparse form.
        D (array): 1D array representing the main diagonal of the degree matrix.
    '''
    hhh, im_brightness = getImage(filename)

    m,n = np.shape(im_brightness)
    im_flat = im_brightness.flatten()
    mn = m*n

    D = np.zeros_like(im_flat)
    W = lil_matrix((mn, mn))
    #print 1

    for i in xrange(len(im_flat)):
        #print i
        indices, distances = getNeighbors(i, radius, m, n)

        areyou = np.exp((-abs(im_flat[indices] - im_flat[i]))/sigma_I - (distances/sigma_d))  #is it squared?
        W[i,indices] = areyou

    D = W.sum(0)
    return sp.csc_matrix(W), D

"""
for i in xrange(mn):
    vecinos = getNeighbors(i, radius, m, n)
    indices = vecinos[0]
    distances = vecinos[1]
    for j in xrange(len(indices)):
        T = indices[j]
        W[i,T] = np.exp(-abs(im_matrix[i] - im_matrix[T])/sigma_I - distances[j]/sigma_d)  #is it squared?
"""

# Problem 4: Implement this function.
def segment(filename="dream.png"):
    '''
    Compute and return the two segments of the image as described in the text.
    Compute L, the laplacian matrix. Then compute D^(-1/2)LD^(-1/2),and find
    the eigenvector corresponding to the second smallest eigenvalue.
    Use this eigenvector to calculate a mask that will be usedto extract
    the segments of the image.
    Inputs:
        filename (string): filename of the image to be segmented
    Returns:
        seg1 (array): an array the same size as img_brightness, but with 0's
                for each pixel not included in the positive
                segment (which corresponds to the positive
                entries of the computed eigenvector)
        seg2 (array): an array the same size as img_brightness, but with 0's
                for each pixel not included in the negative
                segment.
    '''
    hhh, im_brightness = getImage(filename)
    image = np.copy(im_brightness)
    m_,n_ = np.shape(im_brightness)
    W, D = adjacency(filename)
    d_2 = np.power(D,-.5)
    #print "d_2"
    m, n = np.shape(W)
    sp_D = sp.spdiags(D, 0, m, n)
    sp_d_2 = sp.spdiags(d_2, 0, m, n)
    #print "sp_d_2"
    L = sp.csgraph.laplacian(W)

    triple = sp_d_2.dot(L.dot(sp_d_2))
    #print "triple"
    eigvals, eigvecs = eigsh(triple, which = "SM")
    #print "eigvals"
    the_vector = eigvecs[:,1]
    #print the_vector
    negative = eigvecs[:,1]

    the_vector = np.reshape(the_vector,(m_,n_))
    #print "the reshape"
    negative = the_vector.copy()
    the_vector[the_vector>0] = True
    the_vector[the_vector<0] = False
    #print "false"
    negative[negative>0] = False
    negative[negative<0] = True
    #print "true"

    pos = np.multiply(the_vector,image)
    #print "pos"
    neg = np.multiply(negative,image)
    #print "displayPosNeg"
    displayPosNeg(hhh,pos,neg)
    return pos, neg

# Helper function used to convert the image into the correct format.
def getImage(filename='dream.png'):
    '''
    Reads an image and converts the image to a 2-D array of brightness
    values.

    Inputs:
        filename (str): filename of the image to be transformed.
    Returns:
        img_color (array): the image in array form
        img_brightness (array): the image array converted to an array of
            brightness values.
    '''
    img_color = plt.imread(filename)
    img_brightness = (img_color[:,:,0]+img_color[:,:,1]+img_color[:,:,2])/3.0
    return img_color,img_brightness

# Helper function for computing the adjacency matrix of an image
def getNeighbors(index, radius, height, width):
    '''
    Calculate the indices and distances of pixels within radius
    of the pixel at index, where the pixels are in a (height, width) shaped
    array. The returned indices are with respect to the flattened version of the
    array. This is a helper function for adjacency.

    Inputs:
        index (int): denotes the index in the flattened array of the pixel we are
                looking at
        radius (float): radius of the circular region centered at pixel (row, col)
        height, width (int,int): the height and width of the original image, in pixels
    Returns:
        indices (int): a flat array of indices of pixels that are within distance r
                   of the pixel at (row, col)
        distances (int): a flat array giving the respective distances from these
                     pixels to the center pixel.
    '''
    # Find appropriate row, column in unflattened image for flattened index
    row, col = index/width, index%width
    # Cast radius to an int (so we can use arange)
    r = int(radius)
    # Make a square grid of side length 2*r centered at index
    # (This is the sup-norm)
    x = np.arange(max(col - r, 0), min(col + r+1, width))
    y = np.arange(max(row - r, 0), min(row + r+1, height))
    X, Y = np.meshgrid(x, y)
    # Narrows down the desired indices using Euclidean norm
    # (i.e. cutting off corners of square to make circle)
    R = np.sqrt(((X-np.float(col))**2+(Y-np.float(row))**2))
    mask = (R<radius)
    # Return the indices of flattened array and corresponding distances
    return (X[mask] + Y[mask]*width, R[mask])

# Helper function used to display the images.
def displayPosNeg(img_color,pos,neg):
    '''
    Displays the original image along with the positive and negative
    segments of the image.

    Inputs:
        img_color (array): Original image
        pos (array): Positive segment of the original image
        neg (array): Negative segment of the original image
    Returns:
        Plots the original image along with the positive and negative
            segmentations.
    '''
    plt.subplot(131)
    plt.imshow(neg)
    plt.subplot(132)
    plt.imshow(pos)
    plt.subplot(133)
    plt.imshow(img_color)
    plt.show()


if __name__ == "__main__":
    pos, neg= segment()
    print "positive","\n",pos,"\n","negative","\n",neg,
    """
    A = np.array([[0,1,0,0,1,1],[1,0,1,0,1,0],[0,1,0,1,0,0],[0,0,1,0,1,1],[1,1,0,1,0,0],[1,0,0,1,0,0]])
    #print A
    print n_components(A)

    lap = laplacian(A)
    print n_components(lap)
    #print segment()
    #print segment("dreamssX.png")
    #print segment("dream.png")
    "[ 3.93688297  4.73283577  5.67783976 ...,  5.80168915  5.80740309  5.84212637]"
"""
