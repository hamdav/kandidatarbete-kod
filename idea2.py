import time
import numpy as np
import itertools
import sys

from Polynomial import Polynomial, Term
from structConsts import zeropols

if len(sys.argv) > 1:
    n = int(sys.argv[1])
    k = int(sys.argv[2])
else:
    print("Too few arguments")

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))

# choose the smallest number of the first "Box"
# All below this value are 0
# It must be 1 or higher and below 12, because there must be at least three number higher than it left.

def incRecursive(boxes, n, ls, maxBox):
    if n == 0:
        return False
    elif n in ls:
        return incRecursive(boxes, n-1, ls, maxBox-1)
    elif boxes[n] == maxBox:
        boxes[n] = 0
        return incRecursive(boxes, n-1, ls, maxBox)
    else:
        boxes[n] += 1
        return True

def increment(boxes, ls):
    return incRecursive(boxes, len(boxes)-1, ls, len(ls)-1)

def neg(gm):
    r = 2
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

def testIfZeropolsSatisfied(boxes, ls):
    # Test if the zeropols are satisfied
    for key in zeropols:
        #if key[1] != 'VV':
            #continue

        newTerms = dict()

        oldTerms = zeropols[key].getTerms()

        #if boxes == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0]:
            #breakpoint()
        for gms in oldTerms:

            # Switch the Us to the corresponding Us
            # and the Vs 
            newgms = tuple(-ls[boxes[-gm]] if gm < 0 else
                           neg(ls[boxes[neg(gm)]]) for gm in gms)
            sign = +1

            if 0 in newgms or newgms[0] == newgms[1]:
                continue

            elif newgms[1] < newgms[0]:
                newgms = (newgms[1], newgms[0])
                sign = -1

            if newgms in newTerms:
                newTerms[newgms] += sign * oldTerms[gms]
            else:
                newTerms[newgms] = sign * oldTerms[gms]

        allZero = True
        for factor in newTerms.values():
            if factor != 0:
                allZero = False
                break

        if not allZero:
            break # Breaks for k in zeropols
    else: # If for k in zeropols wasn't broken, i.e. if all zeropols were zero
        print("Huzza!")
        print(boxes)

#def testIfSSatisfied(boxes, ls):




iterations = 0
for ls in itertools.combinations(range(1, n+1), k):
    ls = [0, *ls]
    boxes = [0]*(n+1)

    for i, l in enumerate(ls):
        boxes[l] = i


    testIfZeropolsSatisfied(boxes, ls)

    # Yes, this skips the first one but we know that can't work anyway...
    #while increment(boxes, ls):
    while increment(boxes, ls):

        iterations += 1
        # Print every ten million
        if iterations % 1000000 == 0:
            print(iterations)

        testIfZeropolsSatisfied(boxes, ls)




"""
for l1 in range(1, 12):
    for i in range(1, l1):
        boxes[i] = 0
    # Choose the smallest number of the second "Box"
    for l2 in range(l1, 13):
        for zeros in powerset(range(l1+1, l2)):
            for i in range(l1+1, l2):
                boxes[i] = 0 if i in zeros else l1
            # And so on
            for l3 in range(l2, 14):
                for zeros in powerset(range(l2+1, l3)):
                    for ones in powerset([e for e in range(l2+1, l3) if e not in zeros]):
                        for i in range(l2+1, l3):
                            if i in zeros:
                                boxes[i] = 0
                            elif i in ones:
                                boxes[i] = l1
                            else:
                                boxes[i] = l2

                        for l4 in range(l3, 15):
                            for zeros in powerset(range(l3+1, l4)):
                                for ones in powerset([e for e in range(l2+1, l3) if e not in zeros]):
                                    for i in range(l2+1, l3):
                                        if i in zeros:
                                            boxes[i] = 0
                                        elif i in ones:
                                            boxes[i] = l1
                                        else:
                                            boxes[i] = l2
                            # Right. now we choose where the ones that are left go
                            # They can all go into either 0, l1, l2, l3 or l4. 
                            print("HI")
                            """

