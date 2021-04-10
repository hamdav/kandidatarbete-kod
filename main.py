import time
import math
import itertools
import numpy as np
#from quicktions import Fraction

from Polynomial import Polynomial, Term
from structConsts import zeropols, cdict

def neg(gm, n, r):
    if gm < 0:
        if gm < -n + r:
            return gm
        elif gm < (-n + r) // 2:
            return gm + (-n + r) // 2
        else:
            return gm - (-n + r) // 2
    if gm > 0:
        if gm > n - r:
            return gm
        elif gm > (n - r) // 2:
            return gm - (n - r) // 2
        else:
            return gm + (n - r) // 2
    else:
        return 0

def generateS(n, r):
    return Polynomial({(-i, neg(i, n, r)): 1 for i in range(1,n+1)})

def generatePowerOfS(n, r, p, trueFactors=False):
    # If trueFactors, returns the power of S with the actual factors, 
    # which are very large for large p
    # Otherwise returns S scaled so that the factor is 1 for all terms
    termsInS = ((-i, neg(i, n, r)) for i in range(1, n+1))
    return Polynomial([Term(tuple(itertools.chain.from_iterable(gms)), 
                            math.factorial(p) if trueFactors else 1,
                            gmsAreSorted=False)
                       for gms in itertools.combinations(termsInS, p)])

start = time.time()
p = generatePowerOfS(52, 4, 52)
breakpoint()
q = p.reduce3(cdict, zeropols)
print(q)
#s = generateS(15, 3)
#p = s*s
#q = p.reduce3(cdict, zeropols)
#r = q.reduce3(cdict, zeropols)
#print(len(r-q))
end = time.time()
print(f"Time duration: {end - start}")

#for key in zeropols:
    #if key[1] == 'UU':
        #print(zeropols[key])
