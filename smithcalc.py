"""
Qiu Shi Wang
Louis Philippe Ignatieff
Jing Han Sun

Programming final project: Smith normal form user interface
"""

from sympy import *
from smith import *
init_printing(use_unicode=True)

print("Smith normal form calculator")
print("\n")
print("For an integer matrix A, this program will return, in order, a tuple of")
print("matrices D, U and V such that D is in Smith normal form and D=UAV for invertible U and V")
print("\n")
print("Enter an integer matrix in the form of a Python list of its rows, each of which is a Python list")
print("Example input: [[1,3,4],[3,-4,0],[-2,-2,6]]")
while 1 == 1:
    x=input()
    eval("print(smith(Matrix1("+ x +")))")
    print("\n")

#this line allows the console to stay open after the first input()

