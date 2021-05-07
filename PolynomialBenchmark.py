import random
import timeit
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("seaborn-darkgrid")

from Polynomial import Polynomial, Term


def generateRandomPolynomial(n, m, k):
    """
    Returns a random polynomial in a space with
    n grassmann numbers in total,
    m terms, and k factors per term.
    """
    randomTerms = dict()
    for _ in range(m):
        randomTerms[tuple(random.sample(range(1, n+1), k=k))] = random.random()

    return Polynomial(randomTerms)


# Check dependence on m
n = 25
m = 100
k = 10
term_num = 10000
poly_num = 1000
ms = []
termtimes = []
polytimes = []



#while num_iterations > 0:
while m <= 100000 and poly_num > 0 and term_num > 0:
    p=generateRandomPolynomial(n, m, k)
    q=generateRandomPolynomial(n, m, k)
    t = Term((1, 2, 3, 4, 5, 6, 7, 8))
    termtime = timeit.timeit(lambda: p + t, number=term_num)
    termtimes.append(termtime / term_num)
    polytime = timeit.timeit(lambda: q + p, number=poly_num)
    polytimes.append(polytime / poly_num)

    ms.append(m)
    print(f"m: {m}, ttime: {termtime / term_num * 1000:.5f} ms, ptime: {polytime / poly_num * 1000:.5f} ms")

    # get new m
    m = m+500

    # If the times it takes is more than 10 seconds, divide num_iterations by 10
    if termtime > 0.1:
        term_num = term_num // 10
    if polytime > 0.1:
        poly_num = poly_num // 10

#plt.plot(ms, termtimes)
plt.plot(ms, polytimes)
plt.xlabel("Antal termer")
plt.ylabel("Tid (ms)")
plt.title("Tid för att addera två polynom som funktion av antalet termer")
plt.show()
