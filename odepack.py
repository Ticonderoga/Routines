# -*- coding: utf-8 -*-

import numpy as np

#    ________________________
#___/ Matrices tridiagonales \________________________________________

def Thomas_alpha_beta(b,a,c):
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
    alpha = np.ones_like(a)
    beta = np.ones_like(b)
    alpha[0]=a[0]
    n=np.size(a)
    for i in range(1,n) :
        beta[i-1]=b[i-1]/alpha[i-1]
        alpha[i]=a[i]-beta[i-1]*c[i-2]
    
    return alpha,beta
        
def Thomas_y(beta,f) :
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
    y=np.ones_like(f)
    y[0]=f[0]
    n=np.size(f)
    for i in range(1,n) :
        y[i]=f[i]-beta[i-1]*y[i-1]
    return y
    
def Thomas_x(y,alpha,c) :
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
    x=np.ones_like(y)
    x[-1]=y[-1]/alpha[-1]
    n=np.size(y)
    for i in range(n-2,-1,-1) :
        x[i]=(y[i]-c[i]*x[i+1])/alpha[i]
    return x
    
def Solution_Thomas(b,a,c,f) :
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
   
#    _________________
#___/ Méthode d'Euler \________________________________________


#    _____________
#___/ Méthode RK2 \________________________________________

#    _____________
#___/ Méthode RK4 \________________________________________

