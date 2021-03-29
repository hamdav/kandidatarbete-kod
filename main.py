from Polynomial import Polynomial, Term

import time
import numpy as np
#from quicktions import Fraction

from structConsts import zeropols, cdict

def generateS(d):
    return Polynomial({(-i,i): 1 for i in range(1,d+1)})



start = time.time()
s = generateS(14)
p = s*s
q = p.reduce3(cdict, zeropols)
#r = q.reduce3(cdict, zeropols)
#print(len(r-q))
end = time.time()
print(f"Time duration: {end - start}")

for key in zeropols:
    if key[1] == 'UU':
        print(zeropols[key])
