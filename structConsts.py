from Polynomial import Polynomial, Term

import time
import numpy as np

def getZeroPolynomial(c, kind, data):
    # Should be written better!

    n = len(data)
    #breakpoint()
    terms = []

    if kind == 'UU':
        for a in range(n):
            for b in range(a+1, n):
                if round(data[a,b,c-1]) != 0:
                    terms.append(Term(gms=(-(b+1), -(a+1)), factor=(round(data[a,b,c-1]))))
    elif kind == 'UV':
        for a in range(n):
            for b in range(n):
                if round(data[a,b,c-1]) != 0:
                    terms.append(Term(gms=(-(a+1), b+1), factor = (round(data[a,b,c-1]))))
    elif kind == 'VV':
        for a in range(n):
            for b in range(a+1,n):
                if round(data[a,b,c-1]) != 0:
                    terms.append(Term(gms=(a+1, b+1), factor = (round(data[a,b,c-1]))))
    else:
        raise ValueError("Should be one of UU, UV, VV")

    return Polynomial(terms)

def getZeropolsAndCdict(filename):
    """
    Takes a filename containing the structure constants for a group and
    creates the zeropolynomials and cdict for those structure constants
    """

    data = np.load(filename, encoding='bytes')
    n = len(data)

    cdict = {(a,b): [i + 1 for i in np.nonzero(data[a-1,b-1,:])[0]] for a in range(1, n+1) for b in range(1, n+1)}

    zeropols = {(c, kind): getZeroPolynomial(c, kind, data) for c in range(1, n+1) for kind in ['UU', 'UV', 'VV']}

    return zeropols, cdict
