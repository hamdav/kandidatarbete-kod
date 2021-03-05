from Polynomial import Polynomial, Term

import time
import numpy as np
#from quicktions import Fraction

def generateS(d):
    return Polynomial({(-i,i): 1 for i in range(1,d+1)})


data = np.load("StructConstG2.bin", encoding='bytes')
d = len(data)

cdict = {(a,b): [i + 1 for i in np.nonzero(data[a-1,b-1,:])[0]] for a in range(1, d+1) for b in range(1, d+1)}

def getZeroPolynomial(c, kind='UU'):
    # Should be written better!
    #breakpoint()
    terms = []

    if kind == 'UU':
        for a in range(d):
            for b in range(a+1, d):
                if data[a,b,c-1] != 0:
                        terms.append(Term(gms=(-(b+1), -(a+1)), factor=(int(data[a,b,c-1]))))
    elif kind == 'UV':
        for a in range(d):
            for b in range(d):
                if data[a,b,c-1] != 0:
                    terms.append(Term(gms=(-(a+1), b+1), factor = (int(data[a,b,c-1]))))
    elif kind == 'VV':
        for a in range(d):
            for b in range(a+1,d):
                if data[a,b,c-1] != 0:
                    terms.append(Term(gms=(a+1, b+1), factor = (int(data[a,b,c-1]))))
    else:
        raise ValueError("Should be one of UU, UV, VV")

    return Polynomial(terms)

zeropols = {(c, kind): getZeroPolynomial(c, kind=kind) for c in range(1, d+1) for kind in ['UU', 'UV', 'VV']}

start = time.time()
s = generateS(14)
p = s*s*s*s
q = p.reduce3(cdict, zeropols)
#r = q.reduce3(cdict, zeropols)
#print(len(r-q))
end = time.time()
print(f"Time duration: {end - start}")
