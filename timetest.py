# This script takes as command line arguments
# n - An integer and the dimension to test
# p - the power with which to raise S to
#
# It first generates random structure constants with dimension n and
# then checks using the projection method if S^p is equivalent to 0
# under the equivalence relations gotten from the structure constants
# and times how long it takes

import sys
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle

from Polynomial import Polynomial
from Reduction import reduce3
from structConsts import getZeroPolynomial
import fakestructconst

plt.style.use("seaborn-darkgrid")


def generateS(n):
    return Polynomial({(-i, i): 1 for i in range(1,n+1)})

def timeFakeReduction(n, p):

    # Generate the fake structure constants
    data = fakestructconst.GenerateFakeStructconst(n, 0)
    cdict = {(a,b): [i + 1 for i in np.nonzero(data[a-1,b-1,:])[0]] for a in range(1, n+1) for b in range(1, n+1)}
    zeropols = {(c, kind): getZeroPolynomial(c, kind, data) for c in range(1, n+1) for kind in ['UU', 'UV', 'VV']}


    # Start the time
    start = time.time()

    # Generate S
    S = generateS(n)

    # Calculate the power
    Spow = 1
    for _ in range(p):
        Spow = Spow * S

    # Project
    Spow_reduced = reduce3(Spow, cdict, zeropols)

    # End time
    end = time.time()
    duration = end - start
    print("n:", n, "p:", p, "Time taken: ", duration)
    return duration

def timeMulti(N, filename="faketimes.txt"):

    lines = []

    for n in range(6, N):
        line = []
        for p in range(2, n):
            line.append(timeFakeReduction(n, p))
        lines.append(line)

    with open(filename, 'wb') as fp:
        pickle.dump(lines, fp)

def createPlot(filename):

    # Get line data from file
    with open(filename, 'rb') as fp:
        lines = pickle.load(fp)

    # Create a scalar mappable used to create the colorbar and 
    # the colors of the lines on the plot
    cmap = mpl.cm.Oranges
    bounds = [5.5+i for i in range(len(lines)+1)]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='neither')
    scalmap = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

    fig, ax = plt.subplots()

    for line, color in zip(lines, scalmap.to_rgba(range(6, len(lines)+6))):
        #ax.plot(np.linspace(0, 1, len(line)+2)[2:], line)
        ax.plot(range(2, 2+len(line)), line, linewidth=2, color=color)

    ax.set_yscale('log')
    ax.set_xlabel(r'Potens av S: $p$ (1)')
    ax.set_ylabel(r'Körtid: $t$ (s)')
    ax.set_title('Tid för projektion av $S^p$ på slumpade nollpolynom\nmed olika dimensioner')


    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                 ticks=range(6, 6+len(lines)),
                 orientation='vertical',
                 label=r"Dimension: $n$")

    plt.show()


if __name__ == "__main__":

    #timeMulti(12, "12_faketimes.txt")
    createPlot("12_faketimes.txt")
