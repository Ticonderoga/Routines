# -*- coding: iso-8859-15 -*-

from numpy import linspace,meshgrid

epsi=2.
gama=14e-1
Th=1000.
Tc=300.
R=105e-2
nlevel=100
perc=0.01
Th_min=(1.-perc)*Th
Tc_min=(1.+perc)*Tc
DT_min=50.
DT_max=Th_min-Tc_min
teta_max=Th_min/Tc_min
teta_min=Th_min/(Th_min-DT_min)


vec_DT=linspace(DT_min,DT_max,nlevel)
vec_teta=linspace(teta_min,teta_max,nlevel)
X,Y=meshgrid(vec_teta,vec_DT)
Bool1=(Y<Th*(X-1.)/X)
Bool2=Y>Tc*(X-1.)
Bool=Bool1*Bool2
