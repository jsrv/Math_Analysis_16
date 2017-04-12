import numpy as np
from numpy.polynomial import legendre as leg


def GaussQ(f,n):
    return sum(np.multiply(leg.leggauss(n)[1],f(leg.leggauss(n)[0])))


if __name__  == "__main__":

    g = lambda x: np.cos(x)
    h = lambda x: np.abs(x)

    abs_results = []
    cos_results = []

    for i in xrange(10,101,10):
        abs_results.append(GaussQ(h,i))
        cos_results.append(GaussQ(g,i))

    lista = np.array(abs_results) - np.array([1])
    listb =  cos_results -np.array(2*np.sin(1))

    print lista
    print listb

    """
    The reason why is that cosine is easier to aproximate to approximate is because of the polynomial form it has.
    the absolute value has a straight corner that does not allow it.
    """
