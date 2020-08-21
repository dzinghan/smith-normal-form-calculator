"""
Qiu Shi Wang
Louis Philippe Ignatieff
Jing Han Sun

Programming final project: Matrix similarity
"""
from smith import *
from Matrix1 import *

'''Verify if 2 matrices are similar using their Smith Normal form'''

def poly(A, x): 
    '''create the characteristic polynomial and plug in the values'''
    #A = Matrix1(a)
    if A.shape[0] != A.shape[1]: #verify if it's a square matrix
        raise ValueError('matrix is not square')

    I = eye(A.shape[0]) #set up the identity matrix of the right size

    pA = A - x*I

    return pA

def similarity(a, b):
    A = Matrix1(a)
    B = Matrix1(b)
    if A.shape != B.shape: #verify if their sizes match
        raise ValueError('matrices are not the same size')
    elif A.shape[0] != A.shape[1] or B.shape[0] != B.shape[1]:
        raise ValueError('matrix is not square')

    
    valuesA = [] #smith forms of pA
    valuesB = [] #smith forms of pB

    #then compare if valuesA match valuesB

    #if 2 polynomials of degree at most n have the same values on n+1 points,
    #then they are the same polynomial
    
    n = 0 #number of test values must be at least the size of the matrix + 1
    x = -1 #to plus each value of x into A - xI starting from x = 0

    while n <= A.shape[0]:
        x += 1
        #choose values of x such that A - xI is not singular
        #so verify that its rank is equal to its size=
        if poly(A, x).rank() == A.shape[0] and poly(B, x).rank() == B.shape[0]:
            pA = smith(poly(A, x)) #compute the smith form for each given value of x
                                   #between x = 0 to x = n
            valuesA.append(pA[0])

            pB = smith(poly(B, x))
            valuesB.append(pB[0])

            n += 1

    return valuesA == valuesB #if this is true then it's similar


print("Matrix similarity calculator")
print("\n")
print("Enter two square matrices of the same size, A and B, with a comma between them")
print("Example input: [[-1, 6],[-2, 6]],[[3, 0], [0, 2]]")

while 1 == 1:
    x=input()
    tf=bool()
    E=eval("similarity("+x+")")
    if E==True:
        print("The two matrices are similar")
    elif E==False:
        print("The two matrices are not similar")

input()
 #This prevents the console from exiting after the user puts the input

#------------TEST----------
'''
A = [[-1, 6],[-2, 6]]
B = [[3, 0], [0, 2]]

C = [[2, -1], [1, 5]]
D = [[82, 101], [-61, -75]]

E = [[5, 1, 0], [0, 5, 1], [0, 0, 6]]
F = [[5, 0, 0], [0, 5, 1], [0, 0, 6]]


a = similarity(A, B) #they are similar
c = similarity(C, D) #they are similar
e = similarity(E, F) #they are not similar
'''
