"""
Qiu Shi Wang
Louis Philippe Ignatieff
Jing Han Sun

Programming final project: Smith normal form
"""

from sympy import *
from Matrix1 import *

#We want this to return a tuple of three matrices: U, D and V, where D=UAV
def smith(A):
    '''computes the Smith Normal form of a matrix'''
    nrows, ncols = A.shape
    j_t = 0
    leftlist = [] #matrices to left multiply onto A, in order
    rightlist = [] 
    for t in range(nrows):
        
        #step 1: choose pivot
        
        while A.col(j_t) == zeros(nrows, 1) and j_t != ncols-1:
            j_t += 1
        if A.col(j_t) == zeros(nrows, 1) and j_t == ncols-1:
            break
        K = 0

        while A[K, j_t] == 0:
            K += 1
        leftlist.append(exchange_rows(A, t, K))
        A = exchange_rows(A, t, K)*A

        #target row and column
        tr = zeros(1, ncols)
        tr[0, j_t] = A[t, j_t]
        tc = zeros(nrows, 1)
        tc[t,0] = A[t,j_t]
        while A.row(t) != tr or A.col(j_t) != tc:
            
            #step 2: improve pivot
            
            for k in range(t, nrows):
                if A[k,j_t] % A[t,j_t] != 0:
            #We use the Bézout property to find integers s and t such that
            #A[t,j_t] * s + A[k,j_t] * t = b

                    (b,S,T) = egcd(A[t, j_t], A[k, j_t])
                    alpha = A[t, j_t]//b
                    gamma = A[k, j_t]//b

                    L_0 = Matrix1([[S, T], [-gamma, alpha]])
                    Id = eye(nrows)
                    Id[t,t] = S
                    Id[t, k] = T
                    Id[k, t] = -gamma
                    Id[k, k] = alpha

                    leftlist.append(Id)
                    A = Id*A
                
        #step 3: eliminating
                    
            leftlist.append(RRcol(A,  j_t))
            A = RRcol(A, j_t)*A

        
        #step 4: repeat step 2 and 3 for columns
            
            A = A.T() #transpose
            nrows, ncols = A.shape
            t, j_t = j_t, t #transpose
        

            for k in range(t, nrows):
                if A[k, j_t] % A[t, j_t]  != 0:
            #We use the Bézout property to find integers s and t such that
            #A[t, j_t] * s + A[k, j_t] * t = b

                    (b, S, T) = egcd(A[t, j_t], A[k, j_t])
                    alpha = A[t, j_t]//b
                    gamma = A[k, j_t]//b

                    L_0 = Matrix1([[S, T], [-gamma, alpha]])
                    Id = eye(nrows)
                    Id[t,t] = S
                    Id[t, k] = T
                    Id[k, t] = -gamma
                    Id[k, k] = alpha

                    rightlist.append(Id.T()) #append the transpose
                    A = Id*A

            rightlist.append(RRcol(A, j_t).T())
            A = RRcol(A, j_t)*A

            A = A.T() #transpose back
            nrows, ncols = A.shape 
            t, j_t = j_t, t #transpose back

            tr[0, j_t] = A[t, j_t] #check target row and column5
            tc[t, 0] = A[t,j_t]
            
        j_t += 1
        if j_t >= ncols:
            break
        
    #step 5: switch rows of A^T such that the nonzero entries are all on the diagonal

    A = A.T() #transpose
    nrows, ncols = A.shape
    for i in range(nrows):
        if A.row(i) != zeros(1, ncols):
            j = 0
            while A[i, j] == 0:
                j += 1
            rightlist.append(exchange_rows(A, i, j).T())
            A = exchange_rows(A, i, j)*A
    A = A.T()
    nrows, ncols = A.shape
            
    #step 6: ensure that the divisibility criterion is satisfied

    Nc = 0
    while A.col(Nc) != zeros(nrows, 1):
        Nc+=1
        if Nc == ncols:
            break
    Nr = 0
    while A.row(Nr) != zeros(1, ncols):
        Nr+=1
        if Nr == nrows:
            break

    N = min(Nc, Nr) 
    #N is now the total number of nonzero entries/rows/columns
    for i in range(N-1):
        if int(A[i+1, i+1]) % int(A[i, i]) != 0:
            A = A.T() #transpose
            nrows, ncols = A.shape
            Id = eye(nrows)
            Id[i, i+1] = 1
            rightlist.append(Id.T())
            A = Id*A
            beta = gcd(int(A[i, i]), int(A[i+1, i+1]))
            
            A = A.T() #transpose
            nrows, ncols = A.shape
            Id = eye(nrows)
            Id[i, i+1] = int((beta-A[i, i])//(A[i+1, i+1]))
            leftlist.append(Id)
            A = Id*A

            #step 3 again
            leftlist.append(RRcol(A, i))
            A = RRcol(A, i)*A

            A = A.T()
            nrows, ncols = A.shape
            rightlist.append(RRcol(A, i).T())
            A = RRcol(A, i)*A

            A = A.T()
            nrows, ncols = A.shape

    #Step 7: make all the entries in the Smith normal form positive
            #This can be done because the operator of multiplying a row by -1 is invertible
            #over Z. In fact, it is its own inverse.

    Id = eye(nrows)
    for i in range(N):
        if sign(A[i, i]) == -1: #if an entry is negative
            Id[i, i] = -1
    leftlist.append(Id)
    A = Id*A
                 
    L0 = eye(nrows) #Multiply together the (invertible) matrices in the leftlist and the rightlist
                  #to create the matrices U and V
    for L in leftlist:
        L0 = L*L0
    R0 = eye(ncols)
    for R in rightlist:
        R0 = R0*R
    return (A, L0, R0)
                
def RRcol(A, j_t):
    '''Row reduces the j_t-th column, leaving only the uppermost entry nonzero'''
    nrows, ncols  =  A.shape
    i  =  0
    while A[i, j_t] == 0:
        i += 1
    Id = eye(nrows)
    for k in range(i+1, nrows):
        Id[k, i] = -A[k, j_t]//A[i, j_t]
    return Id
    
def exchange_rows(self, i, j):
    '''exchange 2 rows'''
    m, n = self.shape
    Id = eye(m)
    if i == j:
        return Id
    # The elementary matrix which exchanges rows, is the row exchanged idendity
    Id[i, i] = 0
    Id[j, j] = 0
    if i < j:
        i, j = j, i
    Id[i, i-abs(i-j)] = 1
    Id[j, j+abs(i-j)] = 1
    return Id

# Not my code from this line down: found at https://www.techiedelight.com/extended-euclidean-algorithm-implementation/
# This is an implementation of the extended Euclidean algorithm to find the coefficients
# that satisfy Bézout's identity

def egcd(a, b):
    a, b = int(a), int(b)
    if a == 0:
        return (b, 0, 1)
    else:
        gcd, x, y = egcd(b % a, a)
        return (gcd, y - (b//a) * x, x) #The last two numbers are our s and t
