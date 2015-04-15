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
    cdef np.ndarray[FTYPE_t, ndim=1] alpha = np.ones(n, dtype= FTYPE)
    cdef np.ndarray[FTYPE_t, ndim=1] beta = np.ones(n-1, dtype= FTYPE)
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
    cdef np.ndarray[FTYPE_t, ndim=1] y = np.ones(n, dtype= FTYPE)
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
    cdef np.ndarray[FTYPE_t, ndim=1] x = np.ones(n, dtype= FTYPE)
    x[-1]=y[-1]/alpha[-1]
    cdef int i
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
   
