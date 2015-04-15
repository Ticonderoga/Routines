# -*- coding: iso-8859-15 -*-
#      modifié le 22.08.2009 00:39:47
#=================================================================================
#      This program is free software; you can redistribute it and/or modify
#      it under the terms of the GNU General Public License as published by
#      the Free Software Foundation; either version 2 of the License, or
#      (at your option) any later version.
#      
#      This program is distributed in the hope that it will be useful,
#      but WITHOUT ANY WARRANTY; without even the implied warranty of
#      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#      GNU General Public License for more details.
#      
#      You should have received a copy of the GNU General Public License
#      along with this program; if not, write to the Free Software
#      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

import matplotlib.pyplot as plt
import scipy as sc
import numpy as np
from scipy.interpolate import *
from Routines.ListTools import *
from Routines.divers import moypaq,moymob
import matplotlib.font_manager
import itertools
import os,sys
prop = matplotlib.font_manager.FontProperties(size=7)
Ratiofig=21.0/29.7
fig_hauteur=6 #inches
fig_longueur=fig_hauteur/Ratiofig
# on défini la liste des styles qui vont augmenter automatiquement
#
couleurs=['b','g','r','m','c','k']
lignes=[ '-','--', '-.' , ':' ]
markers=['o','v','^','<','>','s','p','*','h','H','+','x','D','d','o','v']
# on combine couleurs et lignes 
lig_style= [x+y for x in lignes for y in couleurs]
colorcycler = itertools.cycle(couleurs).next
linecycler = itertools.cycle(lignes).next
markcycler = itertools.cycle(markers).next
stylecycler= itertools.cycle(lig_style).next
n = itertools.count().next

# cela permet de revenir au début
def resetcycler():
    while stylecycler()<>lignes[-1]+couleurs[-1] :
        pass
#~ 
def merge_mask(*args) :
    """
    merge_mask peremt de fucsionner les masques de plusieurs données
    afin de supprimer les données litigieuses. Soit x,y,T des tableaux de 
    même taille 
    x : un np.array 
    y : un np.ma.array i.e. avec un masque a 
    T : un np.ma.array i.e. avec un masque b
    alors x,y,T=merge_mask(x,y,T) donnera 3 tableaux pour les lesquels on 
    a supprimé les valeurs masquées :
    exemple :
    x=np.arange(10)
    y=np.ma.masked_array(x)
    y[4]=np.ma.masked
    T=3*y+x
    T[2]=np.ma.masked
    a,b,c=merge_mask(x,y,T)
    """
    
    mask_tot=reduce(lambda i,j : i+j,[np.ma.getmaskarray(np.ma.masked_array(el)) for el in args])
    L=[]
    for el in args :        
        if np.ma.is_masked(el) :
            el.mask=mask_tot
        else :
            el=np.ma.masked_array(el,mask=mask_tot,copy=True)
        L.append(el)
    
    return [el.compressed() for el in L]
            
def plotCarto_reg(x,y,T,fignum,Type,levels) :
    x,y,T=merge_mask(x,y,T)
    Tab=np.c_[np.ones_like(x),x,y,x*x,y*y,y*x]
    Coefs=np.linalg.solve(np.dot(Tab.transpose(),Tab),np.dot(Tab.transpose(),T))
    ninterp=100
    xi=np.linspace(x.min(),x.max(),ninterp)
    yi=np.linspace(y.min(),y.max(),ninterp)
    X,Y=np.meshgrid(xi, yi)
    Tcalc=Coefs[0]*np.ones((ninterp,ninterp))+Coefs[1]*X+Coefs[2]*Y+ \
        Coefs[3]*X*X+Coefs[4]*Y*Y+Coefs[5]*X*Y
    plt.figure(fignum,figsize=(fig_longueur,fig_hauteur))
    plt.grid()
    if Type=='Fill' :
        CS1=plt.contourf(X,Y,Tcalc,levels)
    else :
        CS1=plt.contour(X,Y,Tcalc,levels)
        plt.clabel(CS1,fmt = '%2.1f', colors = 'k', fontsize=10)
    vmin=levels.min()
    vmax=levels.max()
    plt.scatter(x,y,s=10,c=T,vmin=vmin,vmax=vmax,cmap=CS1.get_cmap())
    plt.ylabel(r'$\rm{y\ \left[\%\right]}$',fontsize=14)
    plt.xlabel(r'$\rm{x\ \left[\%\right]}$',fontsize=14)
    plt.title(r'$\rm{Temp\'erature\ \left[^\circ C\right]}$',fontsize=14)
    plt.axis([0,1,0,1])
    plt.colorbar()

def plotCarto_spline(x,y,T,fignum,Type,levels) :
    x,y,T=merge_mask(x,y,T)
    ninterp=100
    xi=np.linspace(x.min(),x.max(),ninterp)
    yi=np.linspace(y.min(),y.max(),ninterp)
    X,Y=np.meshgrid(xi, yi)
    #~ [tck1, fp1, ier1, msg1]=bisplrep(x,y,T,kx=2,ky=2,task=0, \
        #~ xb=0,xe=1,yb=0,ye=1,s=0.2,full_output=1)
    [tck1, fp1, ier1, msg1]=bisplrep(x,y,T,kx=2,ky=2,task=0,s=0.2,full_output=1)
    NewT = bisplev(xi,xi,tck1)
    plt.figure(fignum,figsize=(fig_longueur,fig_hauteur))
    plt.grid()
    if Type=='Fill' :
        CS1=plt.contourf(X,Y,NewT.transpose(),levels)
    else :
        CS1=plt.contour(X,Y,NewT.transpose(),levels)
        plt.clabel(CS1,fmt = '%2.1f', colors = 'k', fontsize=10)
    vmin=levels.min()
    vmax=levels.max()
    plt.scatter(x,y,s=10,c=T,vmin=vmin,vmax=vmax,cmap=CS1.get_cmap())
    plt.ylabel(r'$\rm{y\ \left[\%\right]}$',fontsize=14)
    plt.xlabel(r'$\rm{x\ \left[\%\right]}$',fontsize=14)
    plt.title(r'$\rm{Temp\'erature\ \left[^\circ C\right]}$',fontsize=14)
    plt.axis([0,1,0,1])
    plt.colorbar()

def dotprod2(a,b,ax):
    return np.sum(a*b,axis=ax)
    
    
def surfature(X,Y,Z) :
    """surfature -  compute gaussian and mean curvatures of a surface
       [k,h] = surfature(x,y,z), where x,y,z are 2d arrays of points on the
       surface.  k and h are the gaussian and mean curvatures, respectively.
       surfature returns 2 additional arguements,
       [k,h,pmax,pmin] = surfature(...), where pmax and pmin are the minimum
       and maximum curvatures at each point, respectively."""
    # First Derivatives
    [Xv,Xu] = np.gradient(X)
    [Yv,Yu] = np.gradient(Y)
    [Zv,Zu] = np.gradient(Z)
    
    # Second Derivatives
    [Xuv,Xuu] = np.gradient(Xu)
    [Yuv,Yuu] = np.gradient(Yu)
    [Zuv,Zuu] = np.gradient(Zu)
    
    [Xvv,Xuv] = np.gradient(Xv)
    [Yvv,Yuv] = np.gradient(Yv)
    [Zvv,Zuv] = np.gradient(Zv)
    
    # Reshape 2D Arrays into Vectors
    
    Xu = np.array([Xu.flatten('F'),Yu.flatten('F'),Zu.flatten('F')]).transpose()
    Xv = np.array([Xv.flatten('F'),Yv.flatten('F'),Zv.flatten('F')]).transpose()
    Xuu = np.array([Xuu.flatten('F'),Yuu.flatten('F'),Zuu.flatten('F')]).transpose()
    Xuv = np.array([Xuv.flatten('F'),Yuv.flatten('F'),Zuv.flatten('F')]).transpose()
    Xvv = np.array([Xvv.flatten('F'),Yvv.flatten('F'),Zvv.flatten('F')]).transpose()
    
    # First fundamental Coeffecients of the surface (E,F,G
    E = dotprod2(Xu,Xu,1)
    F = dotprod2(Xu,Xv,1)
    G = dotprod2(Xv,Xv,1)
    
    m = np.cross(Xu,Xv,axis=1)
    p = np.sqrt(dotprod2(m,m,1))
    
    n = m/np.transpose(np.tile(p,(3,1)))
    
    # Second fundamental Coeffecients of the surface (L,M,N)
    L = dotprod2(Xuu,n,1)
    M = dotprod2(Xuv,n,1)
    N = dotprod2(Xvv,n,1)
    
    [s,t] = np.shape(Z)
    
    # Gaussian Curvature
    
    K=np.reshape((L*N - M*M)/(E*G - F*F),(s,t),'F')
    
    # Mean Curvature
    H =np.reshape((E*N + G*L - 2*F*M)/(2*(E*G - F*F)),(s,t),'F')
    
    # Principal Curvatures
    Pmax = H + np.sqrt(H*H - K)
    Pmin = H - np.sqrt(H*H - K)
    return [K,H,Pmax,Pmin] 

class RbfEx(object):
    """
    Rbf(*args)

    A class for radial basis function approximation/interpolation of
    n-dimensional scattered data.

    Parameters
    ----------
    *args : arrays
        x, y, z, ..., d, where x, y, z, ... are the coordinates of the nodes
        and d is the array of values at the nodes
    function : str, optional
        The radial basis function, based on the radius, r, given by the norm
        (defult is Euclidean distance); the default is 'multiquadric'::

            'multiquadric': sqrt((r/self.epsilon)**2 + 1)
            'inverse multiquadric': 1.0/sqrt((r/self.epsilon)**2 + 1)
            'gaussian': exp(-(r/self.epsilon)**2)
            'linear': r
            'cubic': r**3
            'quintic': r**5
            'thin-plate': r**2 * log(r)

    epsilon : float, optional
        Adjustable constant for gaussian or multiquadrics functions
        - defaults to approximate average distance between nodes (which is
        a good start).
    smooth : float, optional
        Values greater than zero increase the smoothness of the
        approximation.  0 is for interpolation (default), the function will
        always go through the nodal points in this case.
    norm : callable, optional
        A function that returns the 'distance' between two points, with
        inputs as arrays of positions (x, y, z, ...), and an output as an
        array of distance.  E.g, the default::

            def euclidean_norm(x1, x2):
                return sqrt( ((x1 - x2)**2).sum(axis=0) )

        which is called with x1=x1[ndims,np.newaxis,:] and
        x2=x2[ndims,:,np.newaxis] such that the result is a matrix of the distances
        from each point in x1 to each point in x2.

    Examples
    --------
    >>> rbfi = Rbf(x, y, z, d)  # radial basis function interpolator instance
    >>> di = rbfi(xi, yi, zi)   # interpolated values
    """

    def _euclidean_norm(self, x1, x2):
        return np.sqrt( ((x1 - x2)**2).sum(axis=0) )

    def _function(self, r):
        if self.function.lower() == 'multiquadric':
            return np.sqrt((1.0/self.epsilon*r)**2 + 1)
        elif self.function.lower() == 'inverse multiquadric':
            return 1.0/np.sqrt((1.0/self.epsilon*r)**2 + 1)
        elif self.function.lower() == 'gaussian':
            return exp(-(1.0/self.epsilon*r)**2)
        elif self.function.lower() == 'linear':
            return r
        elif self.function.lower() == 'cubic':
            return r**3
        elif self.function.lower() == 'quintic':
            return r**5
        elif self.function.lower() == 'thin-plate':
            result = r**2 * np.log(r)
            result[r == 0] = 0 # the spline is zero at zero
            return result
        else:
            raise ValueError, 'Invalid basis function name'

    def __init__(self, *args, **kwargs):
        self.xi = np.asarray([np.asarray(a, dtype=np.float_).flatten()
                           for a in args[:-1]])
        self.N = self.xi.shape[-1]
        self.di = np.asarray(args[-1]).flatten()

        assert [x.size==self.di.size for x in self.xi], \
               'All arrays must be equal length'

        self.norm = kwargs.pop('norm', self._euclidean_norm)
        r = self._call_norm(self.xi, self.xi)
        self.epsilon = kwargs.pop('epsilon', r.mean())
        self.function = kwargs.pop('function', 'thin-plate')
        self.smooth = kwargs.pop('smooth', 0.0)
        
        self.A = self._function(r) - np.eye(self.N)*self.smooth
        self.B = np.r_[np.ones((1,self.N)),self.xi]
        self.M = np.r_[np.c_[self.A,self.B.transpose()],np.c_[self.B,np.zeros((3,3))]]
        self.sol = np.linalg.solve(self.M, np.r_[self.di,np.zeros((3,))])
        self.Pcoefs = self.sol[-3:]
        self.nodes = self.sol[:-3]

    def _call_norm(self, x1, x2):
        if len(x1.shape) == 1:
            x1 = x1[np.newaxis, :]
        if len(x2.shape) == 1:
            x2 = x2[np.newaxis, :]
        x1 = x1[..., :, np.newaxis]
        x2 = x2[..., np.newaxis, :]
        return self.norm(x1, x2)
    
    def _CalcAB(self, *args):
        args = [np.asarray(x) for x in args]
        assert all([x.shape == y.shape \
                    for x in args \
                    for y in args]), 'Array lengths must be equal'
        self.xa = np.asarray([a.flatten() for a in args], dtype=np.float_)
        r = self._call_norm(self.xa, self.xi)
        self.LastA = self._function(r)
        self.LastB=np.r_[[np.ones_like(self.xa[1,:])],self.xa]
        return self.LastA,self.LastB

    def __call__(self, *args):
        shp = args[0].shape
        A,B=self._CalcAB(*args)
        return np.dot(np.c_[A,B.transpose()], self.sol).reshape(shp)
    
    def optimize(self,CurvStyle,Tout) :
        from openopt import DFP
        self.N_tour=10
        self.xtour=np.linspace(0,1,self.N_tour)
        self.ytour=self.xtour
        self.Tout=Tout*np.ones_like(self.xtour)
        self.Param_init=self.sol
        self.Atour,self.Btour=self._CalcAB(self.xtour,self.ytour)
        self.N_Curv=10
        self.Xgrid_opt,self.Ygrid_opt=np.meshgrid(np.linspace(0,1,self.N_Curv),np.linspace(0,1,self.N_Curv))
        self.ACurv,self.BCurv=self._CalcAB(self.Xgrid_opt,self.Ygrid_opt)
        Param=self.sol

        def F_opt(Param,X) :
            A,B=self._CalcAB(X[0],X[1])
            
            #print "Tsortie",np.dot(np.c_[A,B.transpose()],Param)
            #print "X",X
            
            return np.dot(np.c_[A,B.transpose()],Param)
            
        def F_cont_tour(Param) :
            return self.Tout-np.dot(np.c_[self.Atour,self.Btour.transpose()],Param)
        
        def F_cont_curv(Param) :
            [K,H,Pmax,Pmin] =surfature(self.Xgrid_opt,self.Ygrid_opt,np.dot(np.c_[self.ACurv,self.BCurv.transpose()],Param).reshape((self.N_Curv,self.N_Curv)))  
            if CurvStyle=="Gaussian" :
                return -K.flatten()
            if CurvStyle=="Mean" :
                return -H.flatten()
        
        def F_cont(Param) :
            return np.r_[F_cont_tour(Param),F_cont_curv(Param)]
            #~ return F_cont_curv(Param)
            
        self.p=DFP(F_opt,self.Param_init,self.xi,self.di,c=F_cont)

def plotCarto_RBF(x,y,T,fignum,Type,levels) :
    x,y,T=merge_mask(x,y,T) 
    vmin=levels.min()
    vmax=levels.max()
    rbf = RbfEx(x, y, T,function='thin-plate')
    ninterp=100
    xi=np.linspace(x.min(),x.max(),ninterp)
    yi=np.linspace(y.min(),y.max(),ninterp)
    X,Y=np.meshgrid(xi, yi)
    Txy=rbf(X,Y)
    plt.figure(fignum)
    CS1=plt.contour(X,Y,Txy,levels)
    plt.clabel(CS1,fmt = '%2.1f', colors = 'k', fontsize=10)
    plt.scatter(x,y,s=10,c=T,vmin=vmin,vmax=vmax,cmap=CS1.get_cmap())
    plt.title(r'$\rm{Temp\'erature\ \left[^\circ C\right]}$',fontsize=14)
    plt.axis([0,1,0,1])
    plt.colorbar()

