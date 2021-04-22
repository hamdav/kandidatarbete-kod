import numbers
import copy
#from quicktions import Fraction
import numpy as np
import scipy.sparse
import scipy.sparse.linalg
from itertools import combinations

def sortWithParity(gms):
    """
    Sorts the iterable gms and keeps track of the parity
    Returns sorted tuple and True if parity of sort is even, False if odd. 

    For now, it uses insertion sort
    """

    parity = True

    l = list(gms)
    sortedIndex = 0     # The list is sorted before this index.

    # When sortedIndex becomes the length of the list, the entire list is sorted.
    while sortedIndex < len(l):

        element = l[sortedIndex]
        i = sortedIndex
        sortedIndex += 1

        while(i>0 and l[i-1]>element):
            l[i]=l[i-1]
            i=i-1
            parity = not parity

        l[i]=element

    return tuple(l), parity



class Term:
    def __init__(self, gms=(), factor=1, gmsAreSorted=True):

        if gmsAreSorted:
            self._grassmanNumbers = gms
            self._factor = factor
        else:
            sortedgms, parity = sortWithParity(gms)
            self._grassmanNumbers = sortedgms
            self._factor = factor if parity else -factor


    def rest(self, skipgms):
        """
        Returns a new term with the grassmannumbers not in skipgms.
        Assumes the removed terms are placed in front and sign is changed thereafter

        Example: (U1 U2 U3).rest((U2)) = -U1 U3
        """
        newGMs = []
        newFactor = self._factor
        for i, gm in enumerate(self._grassmanNumbers):
            if gm in skipgms:
                if i % 2 == 1:
                    newFactor *= -1
            else:
                newGMs.append(gm)

        return Term(gms=newGMs, factor=newFactor)

    def getGMs(self):
        return self._grassmanNumbers

    def getFactor(self):
        return self._factor

    def __mul__(self, other):
        """ 
        Returns self multiplied by other
        """

        if isinstance(other, Term):
            i = ii = 0

            # Parity = True means switch signs, False means don't
            parity = False

            newGms = []

            while (i != len(self._grassmanNumbers) and
                   ii != len(other._grassmanNumbers)):

                a = self._grassmanNumbers[i]
                b = other._grassmanNumbers[ii]

                if a == b:
                    return Term(factor=0)

                if b > a:
                    newGms.append(a)
                    i += 1

                if a > b:
                    newGms.append(b)
                    ii += 1
                    # If b is moved an odd number of steps, change the parity
                    if (len(self._grassmanNumbers) - i) % 2:
                        parity = not parity

            if i != len(self._grassmanNumbers):
                newGms += self._grassmanNumbers[i:]
            elif ii != len(other._grassmanNumbers):
                newGms += other._grassmanNumbers[ii:]

            newFactor = self._factor * other._factor * (-1 if parity else 1)

            return Term(gms=tuple(newGms), factor=newFactor)

        elif isinstance(other, numbers.Number):
            return Term(gms=self._grassmanNumbers, factor=other*self._factor)
        else:
            return NotImplemented


    def __add__(self, other):
        if isinstance(other, Term):
            if self._grassmanNumbers == other._grassmanNumbers:
                return Term(gms=self._grassmanNumbers, factor=self._factor + other._factor)
            else:
                return Polynomial([self, other])
        else:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Term):
            if self._grassmanNumbers == other._grassmanNumbers:
                return Term(gms=self._grassmanNumbers, factor=self._factor - other._factor)
            else:
                return Polynomial([self, -other])
        else:
            return NotImplemented

    def __neg__(self):
        return Term(gms=self._grassmanNumbers, factor=-self._factor)

    def __eq__(self, other):
        """
        Terms that are the same up to a constant numerical factor are equal
        """
        if isinstance(other, Term):
            return self._grassmanNumbers == other._grassmanNumbers
        else:
            return NotImplemented

    # To enable using Term as a key in a dictionary
    def __hash__(self):
        """
        Note again that terms that are the same but with different numerical factors
        are given the same hash
        """
        return hash(self._grassmanNumbers)

    def __repr__(self):
        rv = "+" if self._factor >= 0 else ""
        rv += f"{self._factor}"
        for gm in self._grassmanNumbers:
            if gm < 0:
                rv += f" U{abs(gm)}"
            else:
                rv += f" V{gm}"

        return rv



class Polynomial:

    def __init__(self, terms=[]):
        # self._terms is a dict with tuples of grasmannumbers (i.e. pos or neg ints)
        # as keys and the numerical factors as values

        # The constructor can either take a dictionary formatted as above or
        # a list of Terms

        self._normsquared = None

        if isinstance(terms, dict):
            self._terms = terms

        elif isinstance(terms, list):
            self._terms = dict()

            for term in terms:
                gms = term.getGMs()
                if gms in self._terms:
                    self._terms[gms] += term.getFactor()
                else:
                    self._terms[gms] = term.getFactor()

        else:
            raise ValueError("Polynomial constructor takes only dicts or lists of terms. ")

    def getTerms(self):
        return self._terms

    def clearZeros(self, tol = 1e-7):
        
        """
        Remove all terms with a factor smaller than tol=1e-10
        """
        print("hi")
        delterms = []
        for gmtuple, factor in self._terms.items():
            if abs(factor) < tol:
                delterms.append(gmtuple)

        for gmtuple in delterms:
            del self._terms[gmtuple]

    def scalarProd(self, other):
        """
        Take the scalar product of self with other.
        Scalar product is defined as regular scalar product where the 
        elements of the vector are the factors in front of the terms. 
        For example
        (2 X + 3 XY + 5 Z) . (3 XZ + 10 XY + 2 Z + 5 Y) = 3 * 10 + 5 * 2 = 40
        """

        if len(self) < len(other):
            a = self
            b = other
        else:
            a = other
            b = self

        rv = 0
        for gmtuple, factor in a._terms.items():
            if gmtuple in b._terms:
                rv += factor * b._terms[gmtuple]

        return rv

    def normsquared(self):
        """
        Return the scalar product of self with itself
        """
        if self._normsquared is None:
            self._normsquared = 0
            for _, factor in self._terms.items():
                self._normsquared += factor**2

        return self._normsquared

    def modOut(self, basisVectors):
        """
        Returns polynomial that is orthogonal all of the basis polynomials
        basisVectors MUST be orthogonal
        """
        rv = copy.deepcopy(self)
        for p in basisVectors:
            rv = rv - self.scalarProd(p) / p.normsquared() * p
            
        return rv
            
    def getOrthogonalBasisVectors(self, cdict, zeropols):
        """
        Returns a list of vectors (polynomials) that are mutually orthogonal
        and span the projection of the polynomial self on to the nullspace
        i.e. the space spanned by the polynomials equivalent to 0.
        """

        basisVectors = []
        for gmtuple in self._terms:
            for (a, b) in combinations(gmtuple, 2):
                #breakpoint()
                for c in cdict[abs(a), abs(b)]:
                    if a < 0 and b < 0:
                        kind = 'UU'
                    elif a < 0 and b > 0:
                        kind = 'UV'
                    elif a > 0 and b > 0:
                        kind = 'VV'

                    #basisVector = getZeroPolynomial(c, kind=kind) * term.rest((a, b))
                    #basisVector = basisVector.modOut(basisVectors)
                    restterm = Term(gms=gmtuple, factor=1).rest((a,b))
                    basisVector = (zeropols[(c, kind)] * restterm).modOut(basisVectors)
                    if len(basisVector) != 0:
                        basisVectors.append(basisVector)
        return basisVectors

    def getAllBasisVectors(self, cdict, zeropols):

        """
        Returns a list of all possible basis vectors.
        Does this by simply getting all binomials equivalent to 0 and
        multiplying them with all possible continuations. 

        """

        # NOTE THAT THERE MUST BE AN EQUAL NUMER OF US AND VS IN ALL BASIS
        # POLYNOMIALS!!!!

        basisVectors = set()
        

        # This is ugly, do something better to find the number of factors in each term
        for gmtuple in self._terms:
            m=len(gmtuple)
            break

        # TODO Set in a better way
        n = len(zeropols) // 3

        # n is even.
        for (c, kind), pol in zeropols.items():
            if kind == 'UU':
                for us in combinations(range(-n, 0), m//2 - 2):
                    for vs in combinations(range(1, n+1), m//2):
                        restterm = Term(gms=us+vs, factor=1)
                        basisVector = pol * restterm
                        if len(basisVector) != 0:
                            basisVectors.add(basisVector)
            if kind == 'UV':
                for us in combinations(range(-n, 0), m//2 - 1):
                    for vs in combinations(range(1, n+1), m//2 - 1):
                        restterm = Term(gms=us+vs, factor=1)
                        basisVector = pol * restterm
                        if len(basisVector) != 0:
                            basisVectors.add(basisVector)
            if kind == 'VV':
                for us in combinations(range(-n, 0), m//2):
                    for vs in combinations(range(1, n+1), m//2 - 2):
                        restterm = Term(gms=us+vs, factor=1)
                        basisVector = pol * restterm
                        if len(basisVector) != 0:
                            basisVectors.add(basisVector)

        return basisVectors

    def getBasisVectors(self, cdict, zeropols):

        """
        Returns a list of vectors (polynomials) that are NOT mutually orthogonal
        and span the projection of the polynomial self on to the nullspace
        i.e. the space spanned by the polynomials equivalent to 0.

        Unfortunately, it gives to few vectors so far.
        """

        basisVectors = set()
        for gmtuple in self._terms:
            for (a, b) in combinations(gmtuple, 2):
                #breakpoint()
                for c in cdict[abs(a), abs(b)]:
                    if a < 0 and b < 0:
                        kind = 'UU'
                    elif a < 0 and b > 0:
                        kind = 'UV'
                    elif a > 0 and b > 0:
                        kind = 'VV'

                    restterm = Term(gms=gmtuple, factor=1).rest((a,b))
                    basisVector = zeropols[(c, kind)] * restterm
                    if len(basisVector) != 0:
                        basisVectors.add(basisVector)

        return basisVectors

    def reduce(self, cdict, zeropols):
        bvs = self.getOrthogonalBasisVectors(cdict, zeropols)
        return self.modOut(bvs)

    def reduce2(self, cdict, zeropols):
        breakpoint()
        # Get basisvectors:
        bvs = self.getAllBasisVectors(cdict, zeropols)

        # Convert these basis vectors, 
        # which at the moment are polynomials
        # to numerical sparse scipy vectors (or matrix?)

        # Find all gmtuples of all basisvectors
        gmtuplesToIndices = dict()
        indicesToGMtuples = dict()
        i = 0
        for bv in bvs:
            for gmtuple in bv._terms:
                if not gmtuple in gmtuplesToIndices:
                    gmtuplesToIndices[gmtuple] = i
                    indicesToGMtuples[i] = gmtuple
                    i += 1
        for gmtuple in self._terms:
            if not gmtuple in gmtuplesToIndices:
                gmtuplesToIndices[gmtuple] = i
                indicesToGMtuples[i] = gmtuple
                i += 1

        # Construction of sparse matrix
        row_inds = []
        col_inds = []
        data = []
        for i, bv in enumerate(bvs):
            for gmtuple, factor in bv._terms.items():
                col_inds.append(i)
                row_inds.append(gmtuplesToIndices[gmtuple])
                data.append(factor)

        #breakpoint()
        sparse_mat = scipy.sparse.csc_matrix((data, (row_inds, col_inds)),[len(gmtuplesToIndices), len(bvs)])

        # Construct the solution vector
        sol_vec = np.zeros(len(gmtuplesToIndices))
        for gmtuple, factor in self._terms.items():
            sol_vec[gmtuplesToIndices[gmtuple]] = factor


        # Solve the system of equations
        # TODO: Check if lsmr is faster (see https://stackoverflow.com/questions/12067830/differences-between-scipy-sparse-linalg-lsmr-and-scipy-sparse-linalg-lsqr)
        result = scipy.sparse.linalg.lsqr(sparse_mat, sol_vec)
        x = result[0]
        projected_vector = sparse_mat * x

        #breakpoint()
        # Construct the polynomial
        newTerms = dict()
        for i, factor in enumerate(projected_vector):
            newTerms[indicesToGMtuples[i]] = self._terms.get(indicesToGMtuples[i], 
                                                             0) - factor

        return Polynomial(newTerms)

    def reduce3(self, cdict, zeropols):

        # NOTE THAT THERE MUST BE AN EQUAL NUMER OF US AND VS IN ALL BASIS
        # POLYNOMIALS!!!!

        basisVectors = set()


        # This is ugly, do something better to find the number of factors in each term
        for gmtuple in self._terms:
            n=len(gmtuple)
            break

        # TODO Set in a better way
        dimension=len(zeropols)//3

        gmtuplesToIndices = dict()
        indicesToGMtuples = dict()
        row_index = 0
        col_index = 0

        for gmtuple in self._terms:
            if gmtuple not in gmtuplesToIndices:
                gmtuplesToIndices[gmtuple] = row_index
                indicesToGMtuples[row_index] = gmtuple
                row_index += 1

        row_inds = []
        col_inds = []
        data = []
        # n is even.
        i = 0
        for (c, kind), pol in zeropols.items():
            if kind == 'UU':
                nU = n // 2 - 2
                nV = n // 2
            if kind == 'UV':
                nU = n // 2 - 1
                nV = n // 2 - 1
            if kind == 'VV':
                nU = n // 2
                nV = n // 2 - 2

            for us in combinations(range(-dimension, 0), nU):
                for vs in combinations(range(1, dimension+1), nV):
                    if i % 1000 == 0:
                        print(i)
                    i += 1
                    restterm = Term(gms=us+vs, factor=1)
                    basisVector = pol * restterm
                    if len(basisVector) == 0:
                        continue
                    for gmtuple, factor in basisVector._terms.items():
                        if gmtuple not in gmtuplesToIndices:
                            gmtuplesToIndices[gmtuple] = row_index
                            indicesToGMtuples[row_index] = gmtuple
                            row_index += 1
                        row_inds.append(gmtuplesToIndices[gmtuple])
                        col_inds.append(col_index)
                        data.append(factor)
                    col_index += 1


        sparse_mat = scipy.sparse.csc_matrix((data, (row_inds, col_inds)),[len(gmtuplesToIndices), col_index])

        # Construct the solution vector
        sol_vec = np.zeros(len(gmtuplesToIndices))
        for gmtuple, factor in self._terms.items():
            sol_vec[gmtuplesToIndices[gmtuple]] = factor


        # Solve the system of equations
        # breakpoint()
        result = scipy.sparse.linalg.lsmr(sparse_mat, sol_vec, atol=1e-10, btol=1e-10)
        x = result[0]
        print(f"iteration stop reason: {result[1]}")
        projected_vector = sparse_mat * x

        #breakpoint()
        # Construct the polynomial
        newTerms = dict()
        for i, factor in enumerate(projected_vector):
            newTerms[indicesToGMtuples[i]] = self._terms.get(indicesToGMtuples[i], 
                                                             0) - factor

        return Polynomial(newTerms)
    def __add__(self, other):
        """ 
        Add other and self
        """

        if isinstance(other, Polynomial):
            newTerms = self._terms
            for gmtuple, ofactor in other._terms.items():
                if gmtuple in newTerms:
                    newTerms[gmtuple] += ofactor
                else:
                    newTerms[gmtuple] = ofactor

            rv = Polynomial(newTerms)
            rv.clearZeros()
            return rv
        
        elif isinstance(other, Term):
            newTerms = self._terms.copy()
            ogmtuple = other.getGMs
            if ogmtuple in newTerms:
                newTerms[ogmtuple] += other.getFactor()
            else:
                newTerms[ogmtuple] = other.getFactor()

            rv = Polynomial(newTerms)
            rv.clearZeros()
            return rv
        else:
            return NotImplemented
            
    def __sub__(self, other):
        """ 
        Subtract other and self
        """

        if isinstance(other, Polynomial):
            newTerms = self._terms
            for gmtuple, ofactor in other._terms.items():
                if gmtuple in newTerms:
                    newTerms[gmtuple] -= ofactor
                else:
                    newTerms[gmtuple] = -ofactor

            rv = Polynomial(newTerms)
            rv.clearZeros()
            return rv
        
        elif isinstance(other, Term):
            newTerms = self._terms.copy()
            ogmtuple = other.getGMs
            if ogmtuple in newTerms:
                newTerms[ogmtuple] -= other.getFactor()
            else:
                newTerms[ogmtuple] = -other.getFactor()

            rv = Polynomial(newTerms)
            rv.clearZeros()
            return rv
        
        else:
            return NotImplemented
    def __hash__(self):
        return hash(tuple(sorted(self._terms.items())))

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            newTerms = dict()
            for sgmtuple, sfactor in self._terms.items():
                for ogmtuple, ofactor in other._terms.items():

                    newTerm = Term(gms=sgmtuple, factor=sfactor) * \
                        Term(gms=ogmtuple,factor=ofactor)
                    newGmTuple = newTerm.getGMs()

                    if newGmTuple in newTerms:
                        newFactor = newTerms[newGmTuple] + newTerm.getFactor()
                        del newTerms[newGmTuple]
                    else:
                        newFactor = newTerm.getFactor()

                    if newFactor != 0:
                        newTerms[newGmTuple] = newFactor

            return Polynomial(newTerms)

        elif isinstance(other, Term):

            newTerms = dict()

            for gmtuple, factor in self._terms.items():
                newTerm = Term(gms=gmtuple, factor=factor) * other
                if newTerm.getFactor() != 0:
                    newTerms[newTerm.getGMs()] = newTerm.getFactor()

            return Polynomial(newTerms)

        elif isinstance(other, numbers.Number):

            newTerms = self._terms.copy()

            for gmtuple in newTerms:
                newTerms[gmtuple] *= other

            return Polynomial(newTerms)

        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, Term):

            newTerms = dict()

            for gmtuple, factor in self._terms.items():
                newTerm = other * Term(gms=gmtuple, factor=factor)
                if newTerm.getFactor() != 0:
                    newTerms[newTerm.getGMs()] = newTerm.getFactor()

            return Polynomial(newTerms)

        elif isinstance(other, numbers.Number):
            return self * other

        else:
            return NotImplemented

    def __eq__(self, other):
        return self._terms == other._terms

    def __len__(self):
        return len(self._terms)

    def __repr__(self):
        rv = "P("
        for gmtuple, factor in self._terms.items():
            rv += " +" if factor >= 0 else " "
            rv += f"{factor}"
            for gm in gmtuple:
                if gm < 0:
                    rv += f" U{abs(gm)}"
                else:
                    rv += f" V{gm}"

        rv += " )"
        return rv
