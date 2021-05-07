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

    def _clearZeros(self):
        """
        Removes terms with a factor of 0
        """
        dellist = []
        for gmtuple, factor in self._terms.items():
            if factor == 0:
                dellist.append(gmtuple)

        for gmtuple in dellist:
            del self._terms[gmtuple]

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

    def __neg__(self):
        newTerms = self._terms.copy()

        for gmtuple, factor in newTerms.items():
            newTerms[gmtuple] = -factor

        return Polynomial(newTerms)

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

    def __add__(self, other):
        """ 
        Add other and self
        """

        if isinstance(other, Polynomial):
            if len(self) > len(other):
                newTerms = self._terms.copy()
                oTerms = other._terms
            else:
                newTerms = other._terms.copy()
                oTerms = self._terms

            for gmtuple, ofactor in oTerms.items():
                if gmtuple in newTerms:
                    newTerms[gmtuple] += ofactor
                else:
                    newTerms[gmtuple] = ofactor

            rv = Polynomial(newTerms)
            rv._clearZeros()
            return rv

        elif isinstance(other, Term):
            newTerms = self._terms.copy()
            ogmtuple = other.getGMs()
            if ogmtuple in newTerms:
                newTerms[ogmtuple] += other.getFactor()
            else:
                newTerms[ogmtuple] = other.getFactor()

            rv = Polynomial(newTerms)
            rv._clearZeros()
            return rv
        else:
            return NotImplemented

    def __radd__(self, other):
        return self + other


    def __sub__(self, other):
        """ 
        Return self - other
        """

        return self + (-other)

    def __rsub__(self, other):
        """
        Return other - self
        """

        return (-self) + other

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
