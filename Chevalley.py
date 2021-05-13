# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 16:53:00 2021

@author: soric
"""

import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
"""
To generate structure constants, use GenerateConstants(Type, number).
The type and number can correpsond either to the Cartan classification
(A_l, B_l, ...) or the "physicists" classification (SU(n), SO(n),...).
These are related (for the classical groups) through:
SU(n+1) = A_n
SO(2n) = D_n
SO(2n+1) = B_n
Sp(2n) = C_n

The exceptional groups are G_2, F_4, E_6, E_7 and E_8.

The GenerateStructureConstants function will create numpy-binary files 
with the structure constants, root indices and "c"-constants (i.e 
the constants correpsonding to commutations between e-generators).

If Print is set to true, the program will also save a plot of the c-constants
as a pdf. The appearance of this plot can be changed in the funciton PrintConstants.
"""



def RootFromIndex(Index, FundementalRoots):
    root = np.zeros(len(FundementalRoots[0]))
    for i in range(len(FundementalRoots)):
        root = root + Index[i]*FundementalRoots[i]
    for i in range(len(root)):
        root[i] = int(round(root[i]))
    return root

def isElement(element, lista):
    for i in lista:
        if np.array_equal(i, element):
            return True
    return False

def FindRootIndices(FundementalRoots, r1sqr, r2sqr=None, Type = None):
    dim = len(FundementalRoots)
    I1 = np.zeros(dim, dtype = np.int32)
    I1[0] = 1
    Indexes = np.array([I1])
    for k in range(1,dim):
        I = np.zeros(dim)
        I[k] = 1
        Indexes = np.append(Indexes, [I], axis=0)
    if (Type == 'Sp' or Type == 'C'):
        r2sqr = None
        for i in range(dim - 1):
            I = np.zeros(dim, dtype = np.int32)
            for j in range(i, dim-1):
                I[j] = 2
            I[dim-1] = 1
            Indexes = np.append(Indexes, [I], axis=0)
    i = 0
    while(i < len(Indexes)):
        for j in range(dim):
            Delta = 1
            currentIndex = np.copy(Indexes[i])
            currentIndex[j] = currentIndex[j]  + Delta
            currentRoot = RootFromIndex(currentIndex, FundementalRoots)
            while((np.dot(currentRoot, currentRoot) == r1sqr) or (np.dot(currentRoot, currentRoot) == r2sqr)):
                if not isElement(currentIndex, Indexes):
                    Indexes = np.append(Indexes, [currentIndex], axis=0)
                Delta = Delta + 1
                currentIndex = np.copy(Indexes[i])
                currentIndex[j] = currentIndex[j] + Delta
                currentRoot = RootFromIndex(currentIndex, FundementalRoots)
            
            Delta = -1
            currentIndex = np.copy(Indexes[i])
            currentIndex[j] = currentIndex[j] + Delta
            currentRoot = RootFromIndex(currentIndex, FundementalRoots)
            while((currentIndex[j] >= 0) and ((np.dot(currentRoot, currentRoot) == r1sqr) or (np.dot(currentRoot, currentRoot) == r2sqr))):
                if not isElement(currentIndex, Indexes):
                    Indexes = np.append(Indexes, [currentIndex], axis=0)
                Delta = Delta - 1
                currentIndex = np.copy(Indexes[i])
                currentIndex[j] = currentIndex[j] + Delta
                currentRoot = RootFromIndex(currentIndex, FundementalRoots)
        i = i + 1
    L = len(Indexes)
    for i in range(L):
        Indexes = np.append(Indexes, [-Indexes[i]], axis=0)
                
    return Indexes




def GetSimpleRoots(Type, number):
    if Type == 'G' and number == 2:
        return np.array([[1,0], [-3/2, np.sqrt(3)/2]])
    if Type == 'F' and number == 4:
        A = np.array([2, -2, 0, 0])
        B = np.array([0, 2, -2, 0])
        C = np.array([0, 0, 2, 0])
        D = np.array([-1, -1, -1, 1])
        return np.array([A, B, C, D])
    if Type == 'E':
        E1 = np.array([2, -2, 0, 0, 0, 0, 0, 0])
        E2 = np.array([0, 2, -2, 0, 0, 0, 0, 0])
        E3 = np.array([0, 0, 2, -2, 0, 0, 0, 0])
        E4 = np.array([0, 0, 0, 2, -2, 0, 0, 0])
        E5 = np.array([0, 0, 0, 0, 2, -2, 0, 0])
        E6 = np.array([0, 0, 0, 0, 0, 2, -2, 0])
        E7 = np.array([0, 0, 0, 0, 0, 2, 2, 0])
        E8 = np.array([-1, -1, -1, -1, -1, -1, 1, 1])
        if number == 6:
            return np.array([E3, E4, E5, E6, E7, E8])
        if number == 7:
            return np.array([E2, E3, E4, E5, E6, E7, E8])
        if number == 8:
            return np.array([E1, E2, E3, E4, E5, E6, E7, E8])
    if Type == 'SO':
        if number % 2 == 0:
            Type = 'D'
            number = number//2
        elif number % 2 == 1:
            Type = 'B'
            number = (number-1)//2
    if Type == 'SU':
        Type = 'A'
        number = number - 1
    if Type == 'Sp': 
        if number % 2 == 0:
            Type = 'C'
            number = number//2
        else:
            raise Exception('Error: Sp can only have even numbers')
    l = number
    if Type == 'A':
        SimpleRoots = np.zeros([l, l+1])
        for i in range(l):
            Root = np.zeros(l+1)
            Root[i] = 1
            Root[i+1] = -1
            SimpleRoots[i] = Root
        return SimpleRoots
    if Type == 'B':
        SimpleRoots = np.zeros([l, l])
        for i in range(l-1):
            Root = np.zeros(l)
            Root[i] = 1
            Root[i+1] = -1
            SimpleRoots[i] = Root
        Root = np.zeros(l)
        Root[-1] = 1
        SimpleRoots[-1] = Root
        return SimpleRoots
    if Type == 'C':
        SimpleRoots = np.zeros([l, l])
        for i in range(l-1):
            Root = np.zeros(l)
            Root[i] = 1
            Root[i+1] = -1
            SimpleRoots[i] = Root
        Root = np.zeros(l)
        Root[-1] = 2
        SimpleRoots[-1] = Root
        #print(SimpleRoots)
        return SimpleRoots
    if Type == 'D':
        SimpleRoots = np.zeros([l, l])
        for i in range(l-1):
            Root = np.zeros(l)
            Root[i] = 1
            Root[i+1] = -1
            SimpleRoots[i] = Root
        Root = np.zeros(l)
        Root[-1] = 1
        Root[-2] = 1
        SimpleRoots[-1] = Root
        return SimpleRoots
    raise Exception('Error: Unknown group ' + str(Type) + str(number))
    return -1
    
    



def findElement(element, lista):
    for i in range(len(lista)):
            if np.array_equal(lista[i], element):
                return i
    return None


def Killing(j, i, roots):
    return 4*tij(i,j,roots)/(tij(i,i,roots)*tij(j,j,roots))

def CartanMatrix(i,j, roots):
    return 2*np.dot(roots[i], roots[j])/np.dot(roots[j],roots[j])


def tij(i,j,roots):
    S = 0
    for k in range(len(roots)):
        S = S + CartanMatrix(k, i, roots)*CartanMatrix(k, j, roots)
    return S


def extraSpecialPair(i, roots, fundRoots):
    gamma = roots[i]
    for j in range(len(fundRoots)):
        Index = findElement(gamma - fundRoots[j], roots)
        if Index != None:
            return [j, Index]
    raise Exception("No extraspecial pair found")


def findcab(a, b, roots, fundRoots, heights):
    #a = index of alpha, b = index of beta
    a = int(round(a))
    b = int(round(b))
    #print(str(a) + 'and' + str(b))
    #x = index of xi = alpha + beta
    x = findElement(roots[a] + roots[b], roots)
    
    #is alpha + beta root? If not, constant is zero
    if x == None:
        #print('Not a root!')
        return 0
    
    #Is height(alpha) > height(beta)?
    #If yes, c_alpha,beta = -c_beta,alpha
    if heights[a] > heights[b]:
        return -1*findcab(b, a, roots, fundRoots, heights)
    
    #Is alpha negative?
    #If yes, do 3
    if a >= len(roots)/2:
        
        
        #Is beta also negative?
        #If yes, c_alpha,beta = c_-beta,-alpha
        if b >= len(roots)/2:
            return findcab(b-len(roots)/2, a-len(roots)/2, roots, fundRoots, heights)
        
        
        #Is alpha + beta negative?
        if x >= len(roots)/2:
            return round(Killing(x,x, roots)*findcab(b, x-len(roots)/2, roots, fundRoots, heights)/Killing(a,a, roots))
        
        #If alpha + beta positive:
        return round(Killing(x,x, roots)*findcab(a-len(roots)/2, x, roots, fundRoots, heights)/Killing(b, b, roots))
    
    #If alpha is positive:
    ESP = extraSpecialPair(x, roots, fundRoots)
    #ep = index of epsilon, et = index of eta
    ep = ESP[0]
    et = ESP[1]
    
    #Is epsilon = alpha?
    #If yes, c = r+1 where r greatest integer such that beta - r*alpha is a root (we know that -3<=r<=3)
    if ep == a:
        r = None
        for integer in range(-3, 4):
            if isElement(roots[b] - integer*roots[a], roots):
                r = integer
        if r == None:
            raise Exception("No r found (should be impossible since 0 always works, something is wrong)")
        return r+1
    
    #Is epsilon = beta?
    #If yes, constant = -1
    if ep == b:
        return -1
    
    #Else, a biiig formula
    #r in this case is biggest integer such that eta - r*epsilon is a root (we know that -3<=r<=3)
    r = None
    for integer in range(-3, 4):
        if isElement(roots[et] - integer*roots[ep], roots):
            r = integer
    if r == None:
        raise Exception("No r found (should be impossible since 0 always works, something is wrong)")
    
    
    
    cbep = findcab((ep + len(roots)/2)%len(roots), b, roots, fundRoots, heights)
    caet = findcab((et + len(roots)/2)%len(roots), a, roots, fundRoots, heights)
    caep = findcab((ep + len(roots)/2)%len(roots), a, roots, fundRoots, heights)
    cbet = findcab((et + len(roots)/2)%len(roots), b, roots, fundRoots, heights)
    
    
    """In the big formula, we divide by the killing form of beta - epsilon, 
       and alpha - epsilon. However, it is not certain that these are roots, and if they aren't
       the killing form is undefined. Forutnately, we multiply by c_-epsilon,beta and
       c_-epsilon, alpha. If the roots are undefined, then these are zero. Thus, if we find they are
       zero we do not need to calculate the killing form. """
       
       
    if cbep == 0:
        Kbep = 1
    else:
        bep = findElement(roots[b]-roots[ep], roots)
        Kbep = Killing(bep, bep, roots)
    
    if caep == 0:
        Kaep = 1
    else:
        aep = findElement(roots[a]-roots[ep], roots)
        Kaep = Killing(aep, aep, roots)

    return round((Killing(x,x, roots)/(r+1))*(cbep*caet/Kbep - caep*cbet/Kaep))





def FindCConstants(filename, FundementalRoots, gtype = None):
    r1sqr = None
    r2sqr = None
    for i in range(len(FundementalRoots)):
        rsq = np.dot(FundementalRoots[i], FundementalRoots[i])
        if rsq != r1sqr:
            if r1sqr == None:
                r1sqr = rsq
            elif rsq != r2sqr:
                if r2sqr == None:
                    r2sqr = rsq
        
    Indices = FindRootIndices(FundementalRoots, r1sqr, r2sqr, gtype)
    file = open(filename, 'wb')
    np.save(file, Indices)
    #print(Indices)
    #print(len(Indices))
    roots = np.zeros([len(Indices), len(FundementalRoots[0])])
    heights = np.zeros(len(Indices))
    for i in range(len(Indices)):
        root = np.zeros(len(FundementalRoots[0]))
        for k in range(len(FundementalRoots)):
            root = root + Indices[i,k]*FundementalRoots[k]
        roots[i] = root
        heights[i] = sum(Indices[i])
    k = len(roots)
    #print(roots)
    #print(heights)
    
    constants = np.zeros([k,k])
    for i in range(k):
        #start_time = time.time()
        for j in range(k):
            if constants[j, i] != 0:
                constants[i,j] = -constants[j,i]
            else:
                constants[i,j] = findcab(i, j, roots, FundementalRoots, heights)
            """if constants[j,i] != 0:
                constants[i,j] = -constants[j,i]
            elif constants[ni, nj] != 0:
                constants[i, j] = -constants[ni, nj]
            elif constants[nj, ni] != 0:
                constants[i,j] = constants[nj, ni]
            else:
                constants[i,j] = findcab(i, j, roots, FundementalRoots, heights)"""
        #print("--- %s seconds ---" % (time.time() - start_time))
    #print(constants[0])
    return constants
    

def FindLambdas(roots, FundementalRoots, Indices):
    lambdas = np.array([Indices[0]])
    for i in range(1,len(roots)):
        currentlambda = np.zeros(len(FundementalRoots))
        for j in range(len(FundementalRoots)):
            currentlambda[j] = Indices[i,j]*(np.dot(roots[j],roots[j])/np.dot(roots[i],roots[i]))
        lambdas = np.append(lambdas, [currentlambda], axis=0)
    return lambdas

def FindRoots(FundementalRoots, Indices):
    dim = len(FundementalRoots)
    roots = np.array([FundementalRoots[0]])
    for i in range(1,len(Indices)):
        root = np.zeros(len(FundementalRoots[0]))
        for j in range(dim):
            root = root + Indices[i,j]*FundementalRoots[j]
        roots = np.append(roots, [root], axis=0)
    return roots

def FindStructConstants(roots, constants, lambdas):
    if len(roots) != len(constants):
        print('Säd! rc')
    if len(roots) != len(lambdas):
        print('Säd! rl')
    if len(constants) != len(lambdas):
        print('Säd! cl')
    #print('Hej!')
    dimension = len(lambdas) + len(lambdas[0])
    rank = len(lambdas[0])
    StructureConstants = np.zeros([dimension, dimension, dimension])
    for i in range(dimension):
        for j in range(dimension):
            if (i >= dimension - rank) and (j >= dimension - rank):
                #This means we have two commutators ie zero
                continue
            if (i < dimension - rank) and (j < dimension - rank):
                    #This means we have two non commuters
                    #First we check if i = -j:
                    if (i == int(round(j + (dimension - rank)/2))) or (j == int(round(i + (dimension - rank)/2))):
                        #If yes much weirdness:
                        for k in range(dimension - rank, dimension):
                            StructureConstants[i,j,k] = lambdas[i, (k - dimension + rank)]
                    elif constants[i,j] != 0:
                        #Ie there exists cab, otherwise it's zero
                        if findElement(roots[i] + roots[j], roots) != None:
                                StructureConstants[i,j,findElement(roots[i] + roots[j], roots)] = constants[i,j]
            else:
                #This means we have one h and one e
                #We presume i corresponds to h, otherwise we switch signs
                if i > j:
                    a = i - dimension + rank
                    b = j
                    sign = 1
                else:
                    a = j - dimension + rank
                    b = i
                    sign = -1
                Aab = 2*np.dot(roots[a], roots[b])/np.dot(roots[a], roots[a])
                StructureConstants[i,j,b] = sign*Aab
    return StructureConstants



def GetCsAndPrint(FundementalRoots, groupname, maxCoeff, r1sqr, r2sqr=None):
    Indices = FindRootIndices(FundementalRoots, maxCoeff, r1sqr, r2sqr)
    fileindices = 'Indices' + groupname + '.bin'
    file = open(fileindices, 'wb')
    np.save(file, Indices)
    constants = FindCConstants(fileindices, FundementalRoots)
    
    fileconstants = groupname + 'cab.bin'
    file = open(fileconstants, 'wb')
    np.save(file, constants)
    PrintConstants(constants)
    
def GetSCs(FundementalRoots, groupname):
    indname = 'Indices' + groupname + '.bin'
    constname = groupname + 'cab.bin'
    constants = np.load(constname, encoding='bytes')
    Indices = np.load(indname, encoding = 'bytes')
    roots = FindRoots(FundementalRoots, Indices)
    lambdas = FindLambdas(roots, FundementalRoots, Indices)
    SCs = FindStructConstants(roots, constants, lambdas)
    scname = 'StructConst' + groupname + '.bin'
    file = open(scname, 'wb')
    np.save(file, SCs)
    return SCs

def PrintConstants(constants):
    k = len(constants)


    RorgConstants = np.zeros([k,k])
    for i in range(k):
        if i<k/2:
            ni = int(round((i+k/2)%k))
            #print(i)
            #print(ni)
        else:
            ni = k-i-1
        for j in range(k):
            if j<k/2:
                nj = int(round((j+k/2)%k))
            else:
                nj = k-j-1
            RorgConstants[ni, nj] = np.sign(constants[i,j])*(np.abs(constants[i,j]))**(3/5)
    
    for i in range(k):
        for j in range(k):
            if RorgConstants[i,j] != -RorgConstants[j,i]:
                print(':(')
                print(i)
                print(j)
    
    ax = sb.heatmap(RorgConstants, cmap='PRGn', square=True, cbar=False, xticklabels=False, yticklabels=False)
    ax.invert_xaxis()
    ax.plot()
    plt.savefig('PRGn.pdf')
    
def AllFromScratch(groupname, FundementalRoots, gtype = None):
    Indname = 'Indices' + groupname + '.bin'
    constants = FindCConstants(Indname, FundementalRoots, gtype)
    constname = groupname + 'cab.bin'
    file = open(constname, 'wb')
    np.save(file, constants)
    return GetSCs(FundementalRoots, groupname)

def GenerateStructureConstants(Type, number, Print=False):
    SimpleRoots = GetSimpleRoots(Type, number)
    groupname = Type + str(number)
    SCs = AllFromScratch(groupname, SimpleRoots, Type)
    if Print:
        PrintConstants(np.load(groupname + 'cab.bin', encoding='bytes'))
    return SCs








