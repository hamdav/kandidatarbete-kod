import numpy as np
import scipy.linalg

def generateCoeffMatrix(m, r):
    """
    Returns numpy array like
    1 0 0 ...
    0 0 0 ...
    0 1 0 ...   0
    0 0 0 ...
    0 0 1 ...
    : : :
              1 0 0 0
      0       0 1 0 0
              0 0 1 0
              0 0 0 1
    """
    # r should always be even. 
    assert r % 2 == 0

    # Create a matrix of only zeros with the correct shape
    rv = np.matrix(np.zeros((m, m - r//2)))

    # Put the ones is the correct spots
    for i in range(m - r // 2):
        row_index = 2 * i if i <= r // 2 else r//2 + i
        rv[row_index, i] = 1

    # Return the matrix
    return rv

def generateCoeffMatrix2(T):

    # Remove elements close to zero from T
    T_rounded = np.where(np.abs(T) < 1e-10, 0, T)

    # Find the rank
    r = np.count_nonzero(T_rounded)
    m = T.shape[0]

    # r should always be even. 
    assert r % 2 == 0

    # Create a matrix of only zeros with the correct shape
    rv = np.matrix(np.zeros((m, m - r//2)))

    # Put the ones is the correct spots

    row_index = 0
    for col_index in range(m-r//2):
        rv[row_index, col_index] = 1
        if np.count_nonzero(T_rounded[row_index, :]):
            row_index += 2
        else:
            row_index += 1

    return rv


assert np.array_equal(generateCoeffMatrix(7, 4), np.array(
    [[1, 0, 0, 0, 0],
     [0, 0, 0, 0, 0],
     [0, 1, 0, 0, 0],
     [0, 0, 0, 0, 0],
     [0, 0, 1, 0, 0],
     [0, 0, 0, 1, 0],
     [0, 0, 0, 0, 1]]))
       

# Data has f^c_ab at data[a-1, b-1, c-1], the -1s for zero indexing
data = np.load("StructConstG2.bin", encoding='bytes')

# Round data to int
data = np.rint(data)

# Dimension of lie algebra: n
n = len(data)

# Store the computed matrices
# OBS! SHOULD BE STORED AS np.matrix NOT np.array
coeff_matrices = []
ortho_matrices = []

for c in range(n):

    breakpoint()
    # Find the matrix
    A = np.matrix(data[:, :, c])
    # Ahat = c_2^T O_2 c_1^T O_1 A O_1^T c_1 O_2^T c_2
    for coeff_matrix, ortho_matrix in zip(coeff_matrices, ortho_matrices):
        A = coeff_matrix.T * ortho_matrix * A * ortho_matrix.T * coeff_matrix

    # Find the decomposition
    T, Z = scipy.linalg.schur(A)

    # Number of eigenvalues = rank of matrix
    # Remove elements smaller than 1e-10 in absolute value
    rank = np.count_nonzero(np.where(np.abs(T) < 1e-10, 0, T))
    dim = A.shape[0]

    # Create the coefficient matrix
    #coeff_matrix = generateCoeffMatrix(dim, rank)
    coeff_matrix = generateCoeffMatrix2(T)

    # save the computed matrices
    coeff_matrices.append(coeff_matrix)
    ortho_matrices.append(np.matrix(Z))


