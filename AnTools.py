# -*- coding: utf-8 -*-

import numpy as np
import scipy.signal as sg
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
    alpha = np.ones_like(a.astype('float64'))
    beta = np.ones_like(b.astype('float64'))
    alpha[0]=a[0]
    n=np.size(a)
    for i in range(n-1) :
        beta[i]=b[i]/alpha[i]
        alpha[i+1]=a[i+1]-beta[i]*c[i]
    
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
    y=np.ones_like(f.astype('float64'))
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
    x=np.ones_like(y.astype('float64'))
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

def OdeEuler(f,y0,t):
    """Méthode d'Euler explicite"""
    n=t.size
    dt=t[1]-t[0]
    y=np.empty_like(t)
    y[0]=y0
    for i in xrange(1,n) :
        y[i]=y[i-1]+dt*f(y[i-1],t[i-1])
        
    return y
#    _____________
#___/ Méthode RK2 \________________________________________

def OdeRK2(f,y0,t,Param):
    """Méthode de Runge-Kutta d'ordre 2 :
    le paramétre supplémentaire permet de choisir entre 
    Param = 1 : Méthode d'Euler modifiée
    Param = 0.5 : Méthode du point milieu"""
    
    n=t.size
    dt=t[1]-t[0]
    y=np.empty_like(t)
    y[0]=y0
    if Param == 1 :
        for i in xrange(1,n) :
            dy1=dt*f(y[i-1],t[i-1])
            dy2=dt*f(y[i-1]+dy1,t[i-1]+dt)
            y[i]=y[i-1]+0.5*(dy1+dy2)
        return y
    elif Param == 0.5 :
        for i in xrange(1,n) :
            dy1=dt*f(y[i-1],t[i-1])
            dy2=dt*f(y[i-1]+0.5*dy1,t[i-1]+0.5*dt)
            y[i]=y[i-1]+dy2
        return y
    else :
        print "Paramètre non valide : 0.5 ou 1"    

#    _____________
#___/ Méthode RK4 \________________________________________

def OdeRK4(f,y0,t):
    """
    Méthode de Runge-Kutta d'ordre 4
    """
    n=t.size
    dt=t[1]-t[0]
    y=np.empty_like(t)
    y[0]=y0
    for i in xrange(1,n) :
        dy1=dt*f(y[i-1],t[i-1])
        dy2=dt*f(y[i-1]+dy1/2.,t[i-1]+dt/2.)
        dy3=dt*f(y[i-1]+dy2/2.,t[i-1]+dt/2.)
        dy4=dt*f(y[i-1]+dy3,t[i-1]+dt)
        y[i]=y[i-1]+1./6*(dy1+2.*dy2+2.*dy3+dy4)
    return y

def Tinf_sin(t,Tbase,DT,P) :
    """calcul de Tinf sinus"""
    return Tbase+DT*np.sin(2*np.pi*t/P)

def Tinf_square(t,Tbase,DT,P) :
    """calcul de Tinf carrée"""
    return Tbase+DT*sg.square(2*np.pi*t/P)

def Tinf_saw(t,Tbase,DT,P) :
    """calcul de Tinf dents de scie"""
    return Tbase+DT*sg.sawtooth(2*np.pi*t/P)

def Vin_sin(t,Vo,DV,P) :
    """calcul de Vin sinus"""
    return Vo+DV*np.sin(2*np.pi*t/P)

def Vin_square(t,Vo,DV,P) :
    """calcul de Vin carrée"""
    return Vo+DV*sg.square(2*np.pi*t/P)

def Vin_saw(t,Vo,DV,P) :
    """calcul de Vin dents de scie"""
    return Vo+DV*sg.sawtooth(2*np.pi*t/P)
