'''
Louis Philippe Ignatieff
Qiu Shi Wang
Jing Han Sun
'''

'''Here is our implementation of a Matrix1 class mainly for the purpose of calculating
the Smith normal form. We do not explain all the methods in detail.'''


import copy
class Matrix1(object):
    
    def __init__(self, list):
        '''initialize the matrix as a list of lists, every list, a row'''
        self.list = list
        self.nrows = len(list)
        self.ncols = len(list[0])
        self.shape = (self.nrows, self.ncols)
        
    def is_zero(self):
        '''check for the zero matrix'''
        for row in self:
            for entry in row:
                if entry != 0:
                    return False
        return True
    
    def __setitem__(self,coords,v):
        '''initialize item assignments'''
        i, j = coords
        self.list[i][j] = v
        
    def __getitem__(self,coords):
        '''initialize an index method for matrices'''
        if isinstance(coords, int):
            return self.list[coords]
        else:
            i, j = coords
        return self.list[i][j]

    def row(self, index):
        '''returns the row vector/matrix for a given index'''
        return Matrix1([self.list[index]])

    def T(self):
        '''returns the transpose of the matrix'''
        tr = []
        cols = []
        for i in range(self.ncols):
            for row in self.list:
                cols.append(row[i])
                if len(cols) == self.nrows:
                    tr.append(cols)
                    cols = []
        return Matrix1(tr)
        
    def __repr__(self):
        '''initializes a string representation for matrices'''
        return str(self.list)

    def __str__(self):
        '''initializes the string method using repr'''
        return repr(self)

    def __iter__(self):
        '''initializes an iterable for matrices'''
        return iter(self.list)

    def col(self, index):
        '''returns the column vector/matrix for a given index'''
        return Matrix1([self.T().list[index]]).T()

    def __add__(self, other):
        '''initializes standard matrix addition, entry by entry'''
        if self.shape != other.shape:
            raise ValueError('Matrices of different sizes')
        else:
            NewMatrix = []
            for row_a, row_b in zip(self, other):
                NewRow = []
                for a, b in zip(row_a, row_b):
                    NewRow.append(a+b)
                NewMatrix.append(NewRow)
        return Matrix1(NewMatrix)

    def __mul__(self, other):
        '''initializes standard matrix multiplication and right scalar
        multiplication, here we need to copy the matrix in order
        to keep the original intact, not implemented by us'''
        if type(other) == int: #scalar multiplication to the right
            c = copy.deepcopy(self.list)
            C = Matrix1(c)
            for i in range(self.nrows):
                for j in range(self.ncols):
                    C[i,j] *= other
            return C
        
        elif type(other) == Matrix1 and self.ncols == other.nrows: #matrix multiplication
            P = []
            for i in range(self.nrows):
                row = []
                for j in range(other.ncols):
                    Sum = 0
                    for k in range(self.ncols):
                        Sum += self[i,k]*other[k,j]
                    row.append(Sum)
                P.append(row)
            return Matrix1(P)
        
    def __rmul__(self, other): #scalar multiplication to the left
        '''initializes standard scalar multiplication to the left'''
        if type(other) == int:
            c = copy.deepcopy(self.list)
            C = Matrix1(c)
            for i in range(self.nrows):
                for j in range(self.ncols):
                    C[i,j] *= other
            return C

    def __sub__(self, other):
        '''initializes substraction using scalar multiplication and addition'''
        return self + -1*other

    def __eq__(self, other):
        '''initializes equal sign for matrices (if every entry is equal)'''
        if self.shape != other.shape:
            return False
        else:
            for i in range(self.nrows):
                for j in range(self.ncols):
                    if self[i,j] != other[i,j]:
                        return False
            return True
            

    def __ne__(self, other):
        '''initializes inequal sign for matrices using __eq__'''
        return not self == other

    def rowop1(self, i,j): #All the row operation methods return the matrix
        '''exchange rows elementary operation matrix'''
        #with the row operation done, without modifying self
        E = eye(self.nrows)
        return exchange_rows(E,i,j)*self

    def rowop2(self, i, a):
        '''multiplying row elementary operation matrix'''
        E = eye(self.nrows)
        E[i,i] = a
        return E*self

    def rowop3(self, i, j, a):
        '''subtract rows from each other elementary operation matrix'''
        E = eye(self.nrows)
        E[i,j] = a
        return E*self

    def reduction(self):
        '''reduces to upper triangular using Gauss-Jordan elimination'''
        c = copy.deepcopy(self.list)
        C = Matrix1(c)
        D = zeros(self.nrows, self.ncols)
        
        #Step 1: put all nonzero rows at the top of D.
        #This is equivalent to putting all zero rows at bottom.
        
        k = 0
        for i in range(C.nrows):
            if C.row(i) != zeros(1, C.ncols):
                for j in range(C.ncols):
                    D[k,j] = C[i,j]
            k += 1
            
        #Step 2: find the leftmost nonzero column and a nonzero entry in it: put it on top

        t = 0 #number of iterations, also row of pivot
        j = 0 #column of pivot

        while t != D.nrows-1 and j != D.ncols:
            while j != D.ncols and D.col(j) == 2*D.col(j):
                j += 1
            if j == D.ncols:
                return D
            i = t
            while D[i,j] == 0 and i != D.nrows-1:
                i += 1
            if D[i,j] == 0:
                break
            D = D.rowop1(t,i)
            #Step 3: reduce all rows below it
            for r in range(t+1, D.nrows):
                D = D.rowop3(r, t, -D[r,j]/D[t,j])
            t += 1
            j += 1
        return D
    
    def rank(self):
        '''compute rank of matrix based on row-reduction'''
        R = self.reduction()
        x = 0
        for i in R.list:
            if not Matrix1([i]) == 2*Matrix1([i]):
                x += 1
        return x

def eye(n):
    '''constructs the identity matrix'''
    Id = zeros(n,n)
    for i in range(n):
        Id[i,i] = 1
    return Id
                
def zeros(nrows, ncols):
    '''constructs the zero matrix'''
    M = []
    for i in range(nrows):
        R = []
        for j in range(ncols):
            R.append(0)
        M.append(R)
    return Matrix1(M)

def exchange_rows(self, i, j):
    '''row operation #1'''
    m, n = self.shape
    Id = eye(m)
    if i == j:
        return Id
    # The elementary matrix which exchanges rows, is the row exchanged identity
    Id[i, i] = 0
    Id[j, j] = 0
    if i < j:
        i, j = j, i
    Id[i, i-abs(i-j)] = 1
    Id[j, j+abs(i-j)] = 1
    return Id


#--------TEST---------
'''
M = Matrix1([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
N = Matrix1([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
A = Matrix1([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
H = Matrix1([[1,2],[3,4]])
Z = Matrix1([[0,0],[0,0],[0,0]])
'''
#most important reference : https://igraph.org/python/doc/igraph.datatypes-pysrc.html#Matrix.__ne__


            
