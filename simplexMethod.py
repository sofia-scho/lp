import numpy as np #Life Saver
import sys #For file reader
from math import log10, floor #For rounding, avoid errors from algebraic method

def parseFile():
    """Function to extract data from input file from STDIN"""
    #Initialize variables for matrices
     #B Matrix
    B=[]
    #N Matrix
    N=[]
    #A Matrix
    A=[]
    b=[]
    

    #Read from STDIN and fill matrices
    data = sys.stdin.read()
    data = data.splitlines()
    #Convert each entry to a float type
    dictionary = [[float(num) for num in line.split()] for line in data]

    c=dictionary[0]

    for row in dictionary[1:]:
        i = []
        i += [x for x in row[:-1]]
        A.append(i)
        b.append(row[-1])

    for j, num in enumerate(c):
        N.append(j)
    for k, num in enumerate(A, start = j+1):
        B.append(k)

    return B,N,A,b,c

def sigFigs(num):
    """sigFigs Function to take a number and round it to 7 significant figures"""
    if num==0:
        #ignore
        return 0
    else:
        num = round(num,7-int(floor(log10(abs(num))))-1)
    return num

def primalEntering(zN,B,N):
    """Function to find entering variable for primal simplex"""
    varIndex = 0
    z=[0]*(len(B)+len(N))
    for i in N:
        z[i]=zN[varIndex]
        varIndex+=1
    for i in N:
        if z[int(i)] <0:
            return i

def dualExiting(xB,B,N):
    """Function to find exiting variable for dual simplex"""
    varIndex = 0
    x = [0]*(len(B)+len(N))
    for i in B:
        x[i] = xB[varIndex]
        varIndex+=1
    for i in B:
        if x[int(i)]<0:
            return i

def simplex(A,b,c,B,N):
    """Simplex Function for primal feasible dictionaries"""

    while True:
        AB = A[:,B]
        AN = A[:,N]
        ABInverse = np.linalg.inv(AB)
        xB = np.matmul(ABInverse,b)
        N = [int(num) for num in N]
        B = [int(num) for num in B]
        cN = []
        for i in N:
            if i < len(c):
                cN.append(c[int(i)])
            else:
                cN.append(0)
        cB = []
        for j in B:
            if j < len(c):
                cB.append(c[int(j)])
            else:
                cB.append(0)
        zN = np.matmul(np.transpose(np.matmul(ABInverse,AN)),cB) - cN

        #Check optimality
        if sum(1 for entry in zN if entry >=0) == len(zN):
            objective = np.matmul(np.matmul(np.transpose(cB), ABInverse),b)
            objective = sigFigs(objective)
            #List for coordinates of optimal objective value
            points = []
            for s in range(len(N)):
                if s in B:
                    placeHolder = B.index(s)
                    points.append(xB[placeHolder])
                else:
                    points.append(0)
            print("optimal")
            print(objective)
            for point in points:
                print(point, end = " ")
            return

        #Boo gotta keep pivotting
        entering = primalEntering(zN,B,N)
        ARow = A[:,entering]
        #Choose exiting variable based on entering selection
        deltaXB = np.matmul(ABInverse,ARow)
        deltaXN = 0
        options = []
        for variable, deltaVariable in zip(xB,deltaXB):
            if deltaVariable >0:
                options.append(variable/deltaVariable)
            
        #The LP is unbounded there are no possible pivots
        if not options:
            #Print to STDOUT
            print("unbounded")
            return 

        thisOne = min(options)
        place = 0
        for variable,deltaVariable in zip(xB, deltaXB):
            if deltaVariable > 0 and variable/deltaVariable == thisOne:
                u = place
                break
            place+=1

        exiting = B[u]

        xB = xB - thisOne*deltaXB

        #Pivotting
        B.remove(int(exiting))
        B.append(int(entering))
        N.append(int(exiting))
        N.remove(int(entering))


def dual(A,b,c,B,N,):
    """Dual Simplex Function for feasible dictionaries"""

    while(True):
        N = [int(num) for num in N]
        B = [int(num) for num in B]
        AB = A[:,B]
        AN = A[:,N]
        ABInverse = np.linalg.inv(AB)
        xB = np.matmul(ABInverse,b)
        cN = []
        for i in N:
            if i < len(c):
                cN.append(c[int(i)])
            else:
                cN.append(0)
        cB = []
        for j in B:
            if j < len(c):
                cB.append(c[int(j)])
            else:
                cB.append(0)
        zN = np.matmul(np.transpose(np.matmul(ABInverse,AN)),cB)-cN

        if sum(1 for entry in xB if entry >=0) == len(xB):
            if all(number==0 for number in c):
                return(B,N)
            objective = np.matmul(np.matmul(np.transpose(cB),ABInverse),b)
            objective = sigFigs(objective)
            points = []
            for place in range(len(N)):
                if place in B:
                    here = B.index(place)
                    points.append(xB[here])
                else: 
                    points.append(0)
            print("optimal")
            print(objective)
            for point in points:
                print(point, end = " ")
            return(B,N)

        #Not optimal keep pivoting
        exiting = dualExiting(xB,B,N)
        vector = [0]*len(zB)
        index = B.index(exiting)
        vector[index] = 1
        deltaZN = -np.matmul(np.matmul(np.transpose(AN),np.linalg.inv(np.transpose(AB))),vector)
        #Choose entering given exiting
        options = []
        for variable, deltaVariable in zip(zN, deltaZN):
            if deltaVariable >0:
                options.append(variable/deltaVariable)
        #The LP is infeasible, cannot pivot
        if not options:
            print("infeasible")
            return([],[])

        this = min(options)

        thisOne = 0
        for variable, deltaVariable in zip(zN,deltaZN):
            if deltaVariable > 0 and variable/deltaVariable == this:
                here = thisOne
            thisOne+=1

        entering = N[here]

        zN = zN-this*deltaZN
        variable = this

        #Pivotting
        B.remove(int(exiting))
        B.append(int(entering))
        N.append(int(exiting))
        N.remove(int(entering))
  

def main():
    """Runs either Dual or Primal Simplex Method depending on intial feasibility"""

    #Read the input file and get the matrices
    B,N,A,b,c = parseFile()

    
    #Need identity 
    A=np.array(A)
    rows, columns = A.shape
    identity = np.zeros((rows,rows), float)
    np.fill_diagonal(identity,1)
    A=np.hstack((A, identity))


    #Find if primal or dual feasible and use one of the two functions
    primal = True
    dual = True
    if any(num < 0 for num in b):
        primal = False
    if any(num > 0 for num in c):
        dual = False


    #Primal is feasible
    if primal:
        simplex(A,b,c,B,N)
    elif dual:
        dual(A,b,c,B,N)
    else:
        #Zero it
        temp = [0]*len(N)
        (B,N) = dual(A,b,temp,B,N)
        if B != [] and N != []:
        #Hopefully dual initialized   
            simplex(A,b,c,B,N)


if __name__ == "__main__":
    main()


   
