import numpy as np

data = np.load("StructConstSU5.bin", encoding='bytes')

dim = len(data)

def commutator(A, B):
    return np.matmul(A, B) - np.matmul(B, A)

for a in range(dim):
    for b in range(a+1, dim):
        for c in range(b+1, dim):
            testMat = commutator(data[:,:,a], commutator(data[:,:,b], data[:,:,c])) + \
                      commutator(data[:,:,b], commutator(data[:,:,c], data[:,:,a])) + \
                      commutator(data[:,:,c], commutator(data[:,:,a], data[:,:,b]))
            assert np.all((np.abs(testMat)<1e-10))

"""
for a in range(dim):
    for b in range(dim):
        for c in range(dim):
            for e in range(dim):
                assert 0 == np.dot(data[a,:,e]+data[b,:,e]+data[c,:,e],
                                   data[b,c,:]+data[c,a,:]+data[a,b,:])
"""
