# -*- coding: iso-8859-15 -*-
from numpy import *
from scipy import *
from scipy.sparse import *
from pylab import *

# on cherche à construire une matrice creuse 
# de taille 3000 x 1000 avec des 1, et des -1 à des positions 
# précises

M=40
N=30

Values=zeros((3,max(N,M)))
Values[0,:]=1
Values[1,:]=-1
Values[2,:]=1

Positions=[-10,0,10]
M=spdiags(Values,Positions,M,N)
figure(1)
spy(M.todense())
show()
