# -*- coding: utf-8 -*-

import numpy as np
cimport numpy as np


DTYPE = np.int
ctypedef np.int_t DTYPE_t
FTYPE = np.float
ctypedef np.float_t FTYPE_t

cimport cython
@cython.boundscheck(False)
#    ________________________
#___/ Matrices tridiagonales \________________________________________
# Compilation comme suit 
# cython AnToolsPyx.pyx
# puis 
# gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.7 \
#-I/usr/lib64/python2.7/site-packages/numpy/core/include/ -o AnToolsPyx.so AnToolsPyx.c

def Thomas_alpha_beta(np.ndarray [FTYPE_t, ndim=1] b not None, \
                    np.ndarray [FTYPE_t, ndim=1] a not None, \
                    np.ndarray [FTYPE_t, ndim=1] c not None):
    """ Méthodes de Thomas et Résolution de Ax=f
    f est un vecteur de taille n
    A est composée de 3 diagonales :
        b = diag. inférieure de taille n-1
        a = diag. principale de taille n
        c = diag. supérieure de taille n-1
    
    Calcul des alpha et beta de la méthode de Thomas
    en donnant les 3 diagonales
    exemple :
    alpha,beta=Thomas_alpha_beta(b,a,c):
    """
    cdef int n = a.shape[0]
    cdef np.ndarray[FTYPE_t, ndim=1] alpha = np.empty(n, dtype= FTYPE)
    cdef np.ndarray[FTYPE_t, ndim=1] beta = np.empty(n-1, dtype= FTYPE)
    alpha[0]=a[0]
    cdef int i
    for i in range(n-1) :
        beta[i]=b[i]/alpha[i]
        alpha[i+1]=a[i+1]-beta[i]*c[i]
    
    return alpha,beta
        
def Thomas_y(np.ndarray [FTYPE_t, ndim=1] beta not None, \
                    np.ndarray [FTYPE_t, ndim=1] f not None) :
    """ Méthodes de Thomas et Résolution de Ax=f
    f est un vecteur de taille n
    A est composée de 3 diagonales :
        b = diag. inférieure de taille n-1
        a = diag. principale de taille n
        c = diag. supérieure de taille n-1
    
    Calcul des y selon la méthode de Thomas
    i.e. y solution de Ly=f 
    en donnant beta et f
    exemple :
    y=Thomas_y(beta,f)
    """
    cdef int n = f.shape[0]
    cdef np.ndarray[FTYPE_t, ndim=1] y = np.empty(n, dtype= FTYPE)
    y[0]=f[0]
    cdef int i
    for i in range(1,n) :
        y[i]=f[i]-beta[i-1]*y[i-1]
    return y
    
def Thomas_x(np.ndarray [FTYPE_t, ndim=1] y not None, \
                    np.ndarray [FTYPE_t, ndim=1] alpha not None, \
                    np.ndarray [FTYPE_t, ndim=1] c not None) :
    """ Méthodes de Thomas et Résolution de Ax=f
    f est un vecteur de taille n
    A est composée de 3 diagonales :
        b = diag. inférieure de taille n-1
        a = diag. principale de taille n
        c = diag. supérieure de taille n-1
    
    Calcul des x selon la méthode de Thomas
    i.e. x solution de Ux=y 
    en donnant y,alpha et c
    exemple :
    x=Thomas_x(y,alpha,c)
    """

    cdef int n = y.shape[0]
    cdef np.ndarray[FTYPE_t, ndim=1] x = np.empty(n, dtype= FTYPE)
    x[-1]=y[-1]/alpha[-1]
    cdef int i
    for i in range(n-2,-1,-1) :
        x[i]=(y[i]-c[i]*x[i+1])/alpha[i]
    return x
    
def Solve_bunch_tridiag(np.ndarray [FTYPE_t, ndim=2] tb not None, \
                    np.ndarray [FTYPE_t, ndim=2] ta not None, \
                    np.ndarray [FTYPE_t, ndim=2] tc not None, \
                    np.ndarray [FTYPE_t, ndim=2] tf not None) :
    """ Méthodes de Thomas et Résolution de Ax=f
    f est un vecteur de taille n
    A est composée de 3 diagonales :
        b = diag. inférieure de taille n-1
        a = diag. principale de taille n
        c = diag. supérieure de taille n-1
    
    Calcul des x selon la méthode de Thomas
    exemple :
    x=Solution_Thomas(b,a,c,f)
    """

    cdef int nt = ta.shape[1]
    cdef int n = ta.shape[0]

    cdef np.ndarray[FTYPE_t, ndim=1] alpha = np.empty(nt, dtype= FTYPE)
    cdef np.ndarray[FTYPE_t, ndim=1] beta = np.empty(nt-1, dtype= FTYPE)
    cdef np.ndarray[FTYPE_t, ndim=1] y = np.empty(nt, dtype= FTYPE)
    cdef np.ndarray[FTYPE_t, ndim=1] x = np.empty(nt, dtype= FTYPE)
    cdef np.ndarray[FTYPE_t, ndim=2] tx = np.empty((n,nt), dtype= FTYPE)
    cdef int i
    cdef int idx

    with nogil:
        
        for idx in range(n) :
            alpha[0]=ta[idx,0]
            y[0]=tf[idx,0]
            
    
            for i in range(nt-1) :
                beta[i]=tb[idx,i]/alpha[i]
                alpha[i+1]=ta[idx,i+1]-beta[i]*tc[idx,i]
        
            for i in range(1,nt) :
                y[i]=tf[idx,i]-beta[i-1]*y[i-1]
        
            x[-1]=y[-1]/alpha[-1]
            tx[idx,-1]=x[-1]
            for i in range(nt-2,-1,-1) :
                x[i]=(y[i]-tc[idx,i]*x[i+1])/alpha[i]
                tx[idx,i]=x[i]
            
    return tx



def Solvetridiag(np.ndarray [FTYPE_t, ndim=1] b not None, \
                    np.ndarray [FTYPE_t, ndim=1] a not None, \
                    np.ndarray [FTYPE_t, ndim=1] c not None, \
                    np.ndarray [FTYPE_t, ndim=1] f not None) :
    """ Méthodes de Thomas et Résolution de Ax=f
    f est un vecteur de taille n
    A est composée de 3 diagonales :
        b = diag. inférieure de taille n-1
        a = diag. principale de taille n
        c = diag. supérieure de taille n-1
    
    Calcul des x selon la méthode de Thomas
    exemple :
    x=Solution_Thomas(b,a,c,f)
    """
    cdef int n = a.shape[0]
    cdef np.ndarray[FTYPE_t, ndim=1] alpha = np.empty(n, dtype= FTYPE)
    cdef np.ndarray[FTYPE_t, ndim=1] beta = np.empty(n-1, dtype= FTYPE)
    alpha[0]=a[0]
    cdef np.ndarray[FTYPE_t, ndim=1] y = np.empty(n, dtype= FTYPE)
    y[0]=f[0]
    cdef np.ndarray[FTYPE_t, ndim=1] x = np.empty(n, dtype= FTYPE)
    
    cdef int i
    with nogil:
        for i in range(n-1) :
            beta[i]=b[i]/alpha[i]
            alpha[i+1]=a[i+1]-beta[i]*c[i]
    
        for i in range(1,n) :
            y[i]=f[i]-beta[i-1]*y[i-1]
    
        x[-1]=y[-1]/alpha[-1]

        for i in range(n-2,-1,-1) :
            x[i]=(y[i]-c[i]*x[i+1])/alpha[i]
    
    return x
   
   
def Solution_Thomas(np.ndarray [FTYPE_t, ndim=1] b not None, \
                    np.ndarray [FTYPE_t, ndim=1] a not None, \
                    np.ndarray [FTYPE_t, ndim=1] c not None, \
                    np.ndarray [FTYPE_t, ndim=1] f not None) :
    """ Méthodes de Thomas et Résolution de Ax=f
    f est un vecteur de taille n
    A est composée de 3 diagonales :
        b = diag. inférieure de taille n-1
        a = diag. principale de taille n
        c = diag. supérieure de taille n-1
    
    Calcul des x selon la méthode de Thomas
    exemple :
    x=Solution_Thomas(b,a,c,f)
    """
        
    alpha,beta=Thomas_alpha_beta(b,a,c)
    y=Thomas_y(beta,f)
    x=Thomas_x(y,alpha,c)
    return x
    
    
                
   
def Mean_rings(np.ndarray [FTYPE_t, ndim=2] V not None) :
    cdef int m = V.shape[0]
    cdef int n = V.shape[1]
    cdef int i
    cdef FTYPE_t S1=0.0
    cdef FTYPE_t S2=0.0
    
    for i in range(n) :
        S1=S1+V[1,i]
        S2=S2+V[2,i]
    
    return S1/n,S2/n    


def mul3cols(np.ndarray [FTYPE_t, ndim=2] Bm not None, \
             np.ndarray [FTYPE_t, ndim=2] B not None, \
             np.ndarray [FTYPE_t, ndim=2] Bp not None, \
             np.ndarray [FTYPE_t, ndim=2] T not None) :
    cdef int n = T.shape[0]
    cdef int m = T.shape[1]
    cdef np.ndarray[FTYPE_t, ndim=3] supB = np.array([Bm,B,Bp], dtype= FTYPE)
    cdef np.ndarray[FTYPE_t, ndim=2] result = np.empty((n,m), dtype= FTYPE)
    cdef np.ndarray[FTYPE_t, ndim=2] M = np.empty((n,3), dtype= FTYPE)
#    cdef np.ndarray[FTYPE_t, ndim=2] M2 = np.empty((n,3), dtype= FTYPE)
    
    cdef int j
    cdef int iM
    cdef int jM
    cdef FTYPE_t S=0.0
    for j in range(1,m-1) :
        M=supB[:,:,j].T*T[:,j-1:j+2]
        for iM in range(n) :
            S=0.0
            for jM in range(3) :
                S=S+M[iM,jM]
            result[iM,j]=S
            
    result[:,0]=(supB[1:,:,0].T*T[:,0:2]).sum(-1)
    result[:,-1]=(supB[:-1,:,-1].T*T[:,-2:]).sum(-1)

    return result