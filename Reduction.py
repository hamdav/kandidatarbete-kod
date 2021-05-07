from itertools import combinations
import numpy as np
import scipy.sparse

from Polynomial import Polynomial, Term


def clearZeros(p, tol = 1e-7):

    """
    Remove all terms with a factor smaller than tol
    """
    newTerms = p.getTerms().copy()

    delterms = []
    for gmtuple, factor in newTerms.items():
        if abs(factor) < tol:
            delterms.append(gmtuple)

    for gmtuple in delterms:
        del newTerms[gmtuple]

    return Polynomial(newTerms)


def reduce3(p, cdict, zeropols):

    """
    Takes as arguments:
    * p, a Polynomial. It must have the same number of factors in each term and
      exactly half of them must be Vs, the other half being Us
    * cdict, dictionary that has keys (a, b) and value
      "list of all c such that f^c_ab are nonzero"
    * zeropols, a dictionary that has keys (c, kind) where kind is one of
    'UU', 'UV', 'VV' and values the zeropolynomials f^c_ab X^aX^b according to kind
    """

    # Get the term dictionary of the polynomial
    terms = p.getTerms().copy()

    # Find the number of factors per term
    for gmtuple in terms:
        m=len(gmtuple)
        break

    # Find the dimension of the algebra
    dimension=len(zeropols)//3

    # Initialize some dictionaries.
    # This is the map from rownumber to term and vice versa
    gmtuplesToIndices = dict()
    indicesToGMtuples = dict()
    row_index = 0

    # Go through the terms assigning row numbers to each new grassmanntuple
    for gmtuple in terms:
        if gmtuple not in gmtuplesToIndices:
            gmtuplesToIndices[gmtuple] = row_index
            indicesToGMtuples[row_index] = gmtuple
            row_index += 1


    # Now, construct the matrix of the zeropolynomials
    # This is done as a sparse matrix by specifying the row and column
    # index for each nonzero element.
    row_inds = []
    col_inds = []
    data = []
    col_index = 0

    # Loop through the zeropols and for each one, create a set of
    # new ones by multiplying with all possible grassmannumbers
    for (c, kind), pol in zeropols.items():

        # Depending on the kind of polynomial, a different number of
        # Us and Vs need to be in the last factor.
        if kind == 'UU':
            mU = m // 2 - 2
            mV = m // 2
        elif kind == 'UV':
            mU = m // 2 - 1
            mV = m // 2 - 1
        elif kind == 'VV':
            mU = m // 2
            mV = m // 2 - 2

        # Loop through the possible choices for the last factor
        for us in combinations(range(-dimension, 0), mU):
            for vs in combinations(range(1, dimension+1), mV):
                # Create the last factor
                restterm = Term(gms=us+vs, factor=1, gmsAreSorted=True)
                # Create the basisvector polynomial
                basisVector = pol * restterm
                # If the multiplication resulted in zero,
                # continue to the next pick of Vs
                if len(basisVector) == 0:
                    continue

                # If not, add the row and col indeces and the data to the matrix
                for gmtuple, factor in basisVector._terms.items():
                    if gmtuple not in gmtuplesToIndices:
                        gmtuplesToIndices[gmtuple] = row_index
                        indicesToGMtuples[row_index] = gmtuple
                        row_index += 1
                    row_inds.append(gmtuplesToIndices[gmtuple])
                    col_inds.append(col_index)
                    data.append(factor)
                col_index += 1


    # Construct the sparse matrix
    sparse_mat = scipy.sparse.csc_matrix((data, (row_inds, col_inds)),[len(gmtuplesToIndices), col_index])

    # Construct the solution vector
    sol_vec = np.zeros(len(gmtuplesToIndices))
    for gmtuple, factor in terms.items():
        sol_vec[gmtuplesToIndices[gmtuple]] = factor


    # Solve the system of equations
    result = scipy.sparse.linalg.lsmr(sparse_mat, sol_vec, atol=1e-12, btol=1e-12)
    x = result[0]
    print(f"iteration stop reason: {result[1]}")
    print("1 means x is an approximate solution to Ax = b, i.e. the polynomial is a lin comb of zeropols.\n2 means x approximately solves the least-squares problem, i.e. it is not.")
    projected_vector = sparse_mat * x

    # Construct the polynomial
    newTerms = dict()
    for i, factor in enumerate(projected_vector):
        newTerms[indicesToGMtuples[i]] = terms.get(indicesToGMtuples[i], 0) - factor

    return Polynomial(newTerms)
