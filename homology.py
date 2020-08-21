"""
Qiu Shi Wang
Louis Philippe Ignatieff
Jing Han Sun

Programming final project: Homology of a chain complex
"""

from sympy import *
from smith import *

def Homology(a,b):
    A=Matrix1(a)
    B=Matrix1(b)
    if A.shape[0]!=B.shape[1]:
        return ("Matrices are of the wrong size")
    elif not B*A==zeros(B.shape[0], A.shape[1]):
        return ("BA!=0, not a chain complex")
    k=0
    S=smith(A)[0] #first matrix
    Elemdivisors=[]
    for i in range(min(S.shape[0], S.shape[1])):
        k = S[i][i]
        if k!=0:
            Elemdivisors.append(k)
    answer=["The homology of the complex at the middle term is the direct sum of"]
    if A.shape[0]-A.rank()-B.rank()!=0:
        answer.append(" Z^"+str(A.shape[0]-A.rank()-B.rank()))
    for i in Elemdivisors:
        if i!=1:                     
            answer.append(" Z/"+str(i))
    if len(answer)==1:
        return "The homology of the complex at the middle term is the trivial group."
    ansstring=""
    if len(answer)==2:
        return "The homology of the complex at the middle term is"+answer[1]+"."
    
    for j in range(len(answer)-1):
        ansstring+=answer[j]
    ansstring+=" and"
    ansstring+=answer[len(answer)-1]
    ansstring+="." 

    return ansstring

print("Homology calculator")
print("\n")
print("For a chain complex of finitely generated free abelian groups Z^m -> Z^n -> Z^k where A:Z^m -> Z^n and B:Z^n -> Z^k are homomorphisms represented by matrices, the calculator outputs thehomology group at the middle term H=ker(B)/im(A)")
print("\n")
print("Enter integer matrices A and B with a comma between them. Note that homology is only defined when BA=0.")
print("Example input: [[1,2,5,4],[2,4,10,0],[1,2,5,4]], [[1,0,-1],[0,0,0]]")

while 1 ==1:
    x=input()
    eval("print(Homology("+x+"))")
#This prevents the console from exiting after the user puts the input
