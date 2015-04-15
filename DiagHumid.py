#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Program permettant de tracer un diagramme de Mollier
# pour de l'air ou de l'hydrogène.
# 
# NB :
# air : il faut le Cp de l'air dans Cpa
# H2 : il faut le Cp de H2 dans Cpa
#======================================================

from pylab import draw,figure,plot,isinteractive,ion,ioff,axis,xlabel,ylabel,xticks,yticks,text,title,show,setp,gca
from numpy import array,where,mean,floor,ceil,arange,unique,size,exp,arctan2,pi,select,linspace,logspace,log10
from scipy.optimize import brentq
from string import join
import numpy.ma as Ma
import copy
from proprietes_fluides import *
try: import psyco; psyco.full() 
except: pass

class air :
    """Classe qui permet de définir un air il est caractérisé par 
    * air.car : les 6 caractéristiques [H,na,T,HR,Tr,Th] sous forme d'un 
    tableau masqué T,Tr,Th en °C, na en g eau / g air sec, H en J/kg air sec 
    * air.status : une string
        * indefini
        * defini
        * brouillard
        * saturation
    * air.header : les entêtes i.e. list ["H","na","T","HR","Tr","Th"] 
    ======================================
    Fonctions possibles
    * air.definir permet de définir l'air à l'aide de 2 coordonnées 
    exemple air.definir(["T","na"],[20,0.015])
    le statut passe alors à défini
    * air.ou_non_defini permet de savoir les valeurs definies 
    exemple Tot,index=air.ou_non_defini()
        * Tot le nombre total de car. non définies
        * index les index des valeurs non définies
    * air.ou_defini idem ci-dessus
    exemple Tot,index=air.ou_defini()
    * air.affiche()
    exemple print air.affiche() -> [20000 0.015 -- 0.6 -- 17.5]
        * -- les valeurs masquées
    """
    def __init__(self) :
        self.car=Ma.array([-1.,-1.,-1.,-1.,-1.,-1.],mask=[1,1,1,1,1,1],fill_value=-1e3)
        self.status="indefini"
        self.header=["H","na","T","HR","Tr","Th"]

    def definir(self,headertodef,val) :
        self.car[0:]= Ma.masked
        for i in self.header :
            if headertodef.count(i) > 0 :
                self.car[self.header.index(i)]=float(val[headertodef.index(i)])
            else : pass
        if len(headertodef)== 1 :
            print "Definition impossible car 1 seule coordonee"
            self.status="indefini"
        else :
            self.status="defini"

    def ou_non_defini(self) :
        idx=where(self.car.mask==[1,1,1,1,1,1])
        return [6.0-self.car.count(),idx]

    def ou_defini(self) :
        idx=where(self.car.mask==[0,0,0,0,0,0])
        return [self.car.count(),idx]

    def affiche(self) :
        l=[]
        for i in range(6) :
            if Ma.getmask(self.car)[i] :
                l.append("--")
            else :
                if i==0 :
                    s="%0.f" %(self.car[i])
                elif i==1 :
                    s="%0.3f" %(self.car[i])
                else :
                    s="%0.2f" %(self.car[i])
                l.append(s)
        return "[ "+join(l," , ")+" ]"

def extract(indices,list) :
    """permet d'extraire """
    l = [None]*len(indices)
    for i in range(len(indices)):
        l[i] = list[indices[i]]
    return l
    
def Func_calTR(T,na,HR) : 
    """Fonction de calcul utilisée pour déterminer Tr 
    Func_calTR(T,na,HR) = Pvap(T)*HR-na*Ptot/(RatioM+na)
    """
    return Pvap(T)*HR-na*Ptot/(RatioM+na)

def Func_calTh(T,H,HR) :
    """Fonction de calcul utilisée pour déterminer Th 
    Func_calTh(T,H,HR) = Cpa*T*Ptot+(H+Lv*RatioM)*HR*Pvap(T)+...
    (RatioM*Cpv-Cpa)*T*HR*Pvap(T)-H*Ptot
    """
    return Cpa*T*Ptot+(H+Lv*RatioM)*HR*Pvap(T)+(RatioM*Cpv-Cpa)*T*HR*Pvap(T)-H*Ptot

def Fill_from_Hna(air):
    """Permet de calculer les caractéristiques d'un air à partir des 
    coordonnées de référence H et na le statut de l'air passe alors à 
    brouillard ou saturation si besoin"""
    if air.status=="indefini" :
        print "Air indefini calcul impossible"
    else :
        nondef=air.ou_non_defini() 
        H=air.car[0]
        na=air.car[1]
        T1,T2=Tminopt,air.car[2]
        T3,T4=air.car[4],air.car[2]
        for i in nondef[1][0] :
            if i==2 :
                air.car[i]=(H-na*Lv)/(Cpa+Cpv*na)
                T2,T4=air.car[i],air.car[i]
            elif i==3 : 
                air.car[i]=float(na/(na+RatioM)*Ptot/Pvap(air.car[2]))
                if air.car[i]>1.0 :
                    air.status="brouillard"
                    #print "Air dans une zone de brouillard"
                    T1,T2=air.car[2],Tmaxcrit
            elif i==4 :
                if air.car[3]==1.0 :
                    air.car[i]=air.car[2]
                    air.status="saturation"
                else :
                    air.car[i]=brentq(Func_calTR,T1,T2,args=(air.car[1],1.0))
                    if air.status=="brouillard" : 
                        T3,T4=air.car[2],air.car[i]
                    else :
                        T3,T4=air.car[i],air.car[2]
            elif i==5 :
                if air.car[3]==1.0 :
                    air.car[i]=air.car[2]
                    air.status="saturation"
                else :
                    air.car[i]=brentq(Func_calTh,T3,T4,args=(air.car[0],1.0))
        if air.status=="brouillard" :
            pass
        else :
            air.status="defini"
    return air

def Fill_Air(air):
    """Permet de calculer les caractéristiques d'un air à partir de 2 coordonnées 
    par le jeu des relations on se ramème au cas H et na"""
    estdef=air.ou_defini()
    #print "La Pression est de :",Ptot
    while (estdef[0]<6) :
        Prop_input=estdef[1][0][0:2]
        # on verifie si on est dans une zone possible
        # verif sur H
        if (air.car.mask[0]==0) :
            if (air.car[0]<Hmincrit) :
                air.car[0]=Hmincrit
                Fill_Air(air)
                H=air.car[0]
            else :
                H=air.car[0]
        # verif sur na
        if (air.car.mask[1]==0) :
            if air.car[1]<namincrit :
                air.car[1]=namincrit
                Fill_Air(air)
                na=air.car[1]
            else :
                na=air.car[1]
        # verif sur T
        if (air.car.mask[2]==0) :
            if air.car[2]<Tmincrit :
                air.car[2]=Tmincrit
                Fill_Air(air)
                T=air.car[2]
            else :
                T=air.car[2]
        # pas de verif sur HR car > 1 est possible : brouillard
        if (air.car.mask[3]==0) :   
            HR=air.car[3]
        # verif sur TR
        if (air.car.mask[4]==0) :   
            if air.car[4]<Tmincrit :
                air.car[4]=Tmincrit
                Fill_Air(air)
                TR=air.car[4]
            else :
                TR=air.car[4]
        #verif sur Th
        if (air.car.mask[5]==0) :   
            if air.car[5]<Tmincrit :
                air.car[5]=Tmincrit
                Fill_Air(air)
                Th=air.car[5]
            else :
                Th=air.car[5]
        
        # From H and na
        if (Prop_input==[0,1]).all() :
            Fill_from_Hna(air)
        # From H and T
        elif (Prop_input==[0,2]).all() :
            air.car[1]=float((H-Cpa*T)/(Lv + Cpv*T))
            Fill_from_Hna(air)
        # From H and HR
        elif (Prop_input==[0,3]).all() :
            air.car[2]=brentq(Func_calTh,Tminopt,Tmaxcrit,args=(H,HR))
        # From H and TR
        elif (Prop_input==[0,4]).all() : 
            air.car[1]=float(RatioM*Pvap(TR)/(Ptot-Pvap(TR)))
        # From H and Th
        elif (Prop_input==[0,5]).all() : 
            print "Operation impossible air indefini car H et Th"
            break
        # From na and T
        elif (Prop_input==[1,2]).all() : 
            air.car[0]=float(Cpa*T+na*(Lv+Cpv*T))
        # From na and HR
        elif (Prop_input==[1,3]).all() : 
            air.car[2]=brentq(Func_calTR,Tminopt,Tmaxcrit,args=(na,HR))
        # From na and TR
        elif (Prop_input==[1,4]).all() : 
            print "Operation impossible air indefini car na et TR"
            break
        # From na and Th
        elif (Prop_input==[1,5]).all() : 
            air.car[0] = float(Cpa*Th +RatioM*Pvap(Th)/(Ptot-Pvap(Th))*(Lv + Cpv*Th))
        # From T and HR
        elif (Prop_input==[2,3]).all() : 
            air.car[1]=float(RatioM*Pvap(T)*HR/(Ptot-Pvap(T)*HR))
        # From T and TR
        elif (Prop_input==[2,4]).all() : 
            air.car[1]=float(RatioM*Pvap(TR)/(Ptot-Pvap(TR)))
        # From T and Th
        elif (Prop_input==[2,5]).all() : 
            air.car[0] =float( Cpa*Th +RatioM*Pvap(Th)/(Ptot-Pvap(Th))*(Lv + Cpv*Th))
        # From HR and TR
        elif (Prop_input==[3,4]).all() : 
            air.car[1]=float(RatioM*Pvap(TR)/(Ptot-Pvap(TR)))
        # From HR and Th
        elif (Prop_input==[3,5]).all() : 
            air.car[0] = float(Cpa*Th +RatioM*Pvap(Th)/(Ptot-Pvap(Th))*(Lv + Cpv*Th))
        # From TR and Th
        elif (Prop_input==[4,5]).all() :
            air.car[0] = float(Cpa*Th +RatioM*Pvap(Th)/(Ptot-Pvap(Th))*(Lv + Cpv*Th))
        else : pass
        estdef=air.ou_defini()
    return air

def Return_BBox(air_bg,air_hd) :
    """Fonction qui retourne les caractéristiques des airs aux 
    4 coins du diagramme à partir de l'air en bas à gauche et 
    celui en haut à droite"""
    air_hg=air()
    air_hg.definir(["na","H"],[air_bg.car[1],air_hd.car[0]-(air_hd.car[1]-air_bg.car[1])*Lv])
    Fill_Air(air_hg)
    air_bd=air()    
    air_bd.definir(["na","H"],[air_hd.car[1],air_bg.car[0]-(air_bg.car[1]-air_hd.car[1])*Lv])
    Fill_Air(air_bd)
    a1=copy.deepcopy(air_bg)
    a3=copy.deepcopy(air_hd)
    return a1,air_bd,a3,air_hg

def resple(li,nval,type) :
    """fonction qui permet de resampler li selon un nombre de val 
    nval avec un Delta dv"""
    if type=="H" :
        Vecdv=[1000,5000,10000,15000,20000,25000,30000,50000]
    elif type=="na" :
        Vecdv=[5e-4,1e-3,2e-3,1e-2,0.1,0.5]
    elif type=="T" :
        Vecdv=[0.1,0.5,1,2,5,10,100]
    elif type=="HR" :
        Vecdv=[0.01,0.05,0.1,0.2]
    liinit=li
    i=0 
    while len(li)>=nval :
        #li=unique(liinit//Vecdv[i])*Vecdv[i]
        li=unique(floor(liinit/Vecdv[i])*Vecdv[i])
        i=i+1
    return list(li)

def Return_levels_HnaTHR(ain,aout) :
    """Fonction qui renvoie les niveaux de H,na,T et HR à afficher"""
    a1,a2,a3,a4=Return_BBox(ain,aout)  
    namin=max(floor(a1.car[1]*1000)/1000,namincrit)
    namax=ceil(a2.car[1]*1000)/1000
    #print namin,namax
    #print "A1",a1.affiche()
    #print "A2",a2.affiche()
    #print "A3",a3.affiche()
    #print "A4",a4.affiche()
    
    if floor(a1.car[2])<floor(a2.car[2]) :
        Tmin=floor(a1.car[2])
        a1.definir(["na","T"],[namin,Tmin])
        Fill_Air(a1)
        a2.definir(["na","H"],[namax,a1.car[0]-(namin-namax)*Lv])
        Fill_Air(a2)
        #print "Cas 1 T1<T2 Tmin=",Tmin     
    else :
        Tmin=floor(a2.car[2])
        a2.definir(["na","T"],[namax,Tmin])
        Fill_Air(a2)
        a1.definir(["na","H"],[namin,a2.car[0]-(namax-namin)*Lv])
        Fill_Air(a1)    
        #print "Cas 2 T1>T2 Tmin=",Tmin
        
    if ceil(a3.car[2])<ceil(a4.car[2]) :
        Tmax=ceil(a4.car[2])
        a4.definir(["na","T"],[namin,Tmax])
        Fill_Air(a4)
        a3.definir(["na","H"],[namax,a4.car[0]-(namin-namax)*Lv])
        Fill_Air(a3)
    else :
        Tmax=ceil(a3.car[2])
        a3.definir(["na","T"],[namax,Tmax])
        Fill_Air(a3)
        a4.definir(["na","H"],[namax,a3.car[0]-(namax-namin)*Lv])
        Fill_Air(a4)
    
    HRmin=floor(a4.car[3]*100)/100
    HRmax=ceil(min(1.0,a2.car[3])*100)/100
    Hmin=floor(a1.car[0]*1e3)/1e3
    Hmax=ceil(a3.car[0]*1e3)/1e3    
    na_base=arange(0.0,0.5,0.0005)
    HR_base=arange(0.0,1.001,0.01)
    T_base=arange(Tmincrit,Tmaxcrit+0.01,0.1)
    H_base=arange(-29000,1.9e6+1.0,1e3)
    list_H=resple(H_base[where((H_base>=Hmin) & (H_base<=Hmax))],50,"H")
    #list_na=resple(na_base[where((na_base>=namin) & (na_base<=namax))],100,"na")
    #list_na=arange(namin,namax+0.001,0.001)
    list_na=list(linspace(namin,namax,num=100,endpoint=True))
    #print type(list_na)
    if list_na[0]==0.0 :
        list_na.remove(0)
    list_T=resple(T_base[where((T_base>=Tmin) & (T_base<=Tmax))],30,"T")
    list_HR=resple(HR_base[where((HR_base>=HRmin) & (HR_base<=HRmax))],30,"HR")
    if list_HR[0]==0.0 :
        list_HR.remove(0)   
    #print "A1",a1.affiche()
    #print "A2",a2.affiche()
    #print "A3",a3.affiche()
    #print "A4",a4.affiche()
    return list_H,list_na,list_T,list_HR,[a1,a2,a3,a4]

def Transf_xy(x,y) :
    """Effectue une tranformation d'axe
    de na,H en x,y"""
    y=y-x*Lv
    return x,y

def Return_Hs(lh,lna) :
    """Permet de calculer des isoH avec
    * air1 et air2 en entree
    * n = nombre d'isoH"""
    list_H=[]
    list_na=[]
    for i in lh :
        list_H.append(i)
        list_H.append(i)
        list_na.append(lna[0])
        list_na.append(lna[-1])
    return list_na,list_H

def Return_Ts(lt,lna) :
    """Permet de calculer des isoT avec
    * air1 et air2 en entree
    * n = nombre d'isoT"""
    list_H=[]
    list_na=[]
    Atmp=air()
    for i in lt :
        Atmp.definir(["T","na"],[i,lna[0]])
        Fill_Air(Atmp)
        list_H.append(Atmp.car[0])
        Atmp.definir(["T","na"],[i,lna[-1]])
        Fill_Air(Atmp)
        list_H.append(Atmp.car[0])  
    list_na.append(lna[0])
    list_na.append(lna[-1])
    return list_na,list_H
    
def Return_HRs(lna,lhr) :
    """Permet de calculer des isoHR avec
    * """
    list_H=[]
    metalist_H=[]
    list_na=[]
    Atmp=air()
    namin=lna[0]
    namax=lna[-1]
    A=10.
    e=exp(A)
    C=(namax-namin)/(e-1)
    Vecna=C*(exp(A*(lna-namin)/(namax-namin))-1)+namin
    #Vecna=logspace(log10(namin),log10(namax),num=100,endpoint=True)
    for i in lhr :
        for j in Vecna :
            Atmp.definir(["HR","na"],[float(i),float(j)])
            Fill_Air(Atmp)
            list_H.append(Atmp.car[0])
        metalist_H.append(list_H)
        list_H=[]
    return Vecna,metalist_H
    
def Plotlines(lna,lH,line) :
    """Permet de tracer les lignes sur le diagramme avec 
    * lna liste de taille 2
    * lH liste de taille 2*n avec n le nombre de lignes à tracer
    * line style de la ligne"""
    
    slH=size(lH)
    for i in arange(0,slH,2) :
        nx=[lna[0],lna[-1]]
        ny=[Transf_xy(lna[0],lH[i])[1],Transf_xy(lna[-1],lH[i+1])[1]]
        #nx=array(nx)
        #ny=array(ny)
        #nx=select([nx<x1,nx>x3],[x1,x3], default = nx)
        #ny=select([ny<y1,ny>y3],[y1,y3], default = ny)
        #print nx
        #print ny       
        #plot(nx.tolist(),ny.tolist(),**line)
        plot(nx,ny,**line)

def PlotHR(lna,metalH,lhr,line) :
    """Permet de tracer des isoHRs sur le diagramme avec   
    * lna liste de taille nb
    * metalH liste (n élements) de listes (nb élements)
        n nombre d'isoHR
        nb nombre de points sur une HR
    * line style de la ligne"""
    ft_HR=dict(horizontalalignment='center',
                  verticalalignment='center',color='k',
                  size=8)
    smlH=len(metalH)
    i_coup=(len(lna)-1)-smlH
    for i in arange(0,smlH,1) :
        lres=zip(*map(Transf_xy,lna,metalH[i]))
        plot(lres[0],lres[1],**line)
        if (i+1)%2==0 :
            DX=(lres[1][i_coup+1+i]-lres[1][i_coup-1+i])*Ratiofig
            DY=(lres[0][i_coup+1+i]-lres[0][i_coup-1+i])*RatioAxes
            angle=arctan2(DX,DY)/pi*180
            if (lres[0][i_coup+i]<=x3*0.95) and (lres[0][i_coup+i]>=x1*0.95) \
            and (lres[1][i_coup+i]<=y3*0.95) and (lres[1][i_coup+i]>=y1*0.95) :
                text(lres[0][i_coup+i],lres[1][i_coup+i],"HR="+str(lhr[i])
                      ,rotation=angle,**ft_HR)
            else : pass     

def Init_diag(Ain,Aout,fignum):
    """ Procédure d'initialisation du diagramme et qui permet de fixer 
    les limites du diagramme. Par exemple, on fixe ainsi dans une variable 
    globale la pression totale. Voilà un exemple d'utilisation :
    
    n=1
    A_bg=Diag.air()
    A_bg.definir(["T","na"],[20.0,0.00001])
    A_hd=Diag.air()
    A_hd.definir(["T","na"],[80.0,0.6])
    LH,Lna,LT,LHR,Lair=Diag.Init_diag(A_bg,A_hd,n)
    
    avec A_bg l'air en bas à gauche et A_hd l'air en haut à droite
    le dernier paramètre ici n=1 spécifie le numéro de la figure matplotlib. 
    LH,Lna,LT,LHR,Lair sont les sorties et correpondent aux niveaux sur
    les différentes coordonnées. Lair est une liste de taille 4 qui contient 
    les airs des 4 coins du diagramme. 
    
    si l'on met n=-1 alors on ne fait que des calculs sur des airs humides
    et l'appel se fait ainsi :
    
    Diag.Init_diag(A_bg,A_hd,-1)
    """
    global Ptot,RatioM,Tmincrit,Cpa,Cpv,Lv,namincrit,Hmincrit,Tmaxcrit,Tminopt,Rgaz,MH2O,Mair
    global Ratiofig,fig_hauteur,fig_longueur,RatioAxes,lineprops,lineprops_HR,x1,y1,x3,y3
    Rgaz=8.314472
    Ptot=101325.0
    MO=15.9994
    MN=14.0067
    MAr=39.948
    Mair=0.78084*2*MN+0.20946*2*MO+0.009340*MAr
    MH=1.00794
    MH2O=2*MH+MO
    RatioM=MH2O/Mair
    Tmincrit=-30.0
    Cpa=Cpas(C2K(Tmincrit))
    Cpv=Cpvap(C2K(0.1))
    Lv=ChLv(C2K(0.1))
    namincrit=RatioM*Pvap(Tmincrit)/(Ptot-Pvap(Tmincrit))
    Hmincrit=Cpa*Tmincrit+namincrit*(Lv+Cpv*Tmincrit)
    Tmaxcrit=260.0
    Tminopt=-30.0
    Tinit=20.0
    Cpa=Cpas(C2K(Tinit))
    Cpv=Cpvap(C2K(Tinit))
    Lv=ChLv(C2K(Tinit))
    maskin=Ma.getmask(Ain.car)
    maskout=Ma.getmask(Aout.car)
    
# Il s'agit ici de fixer les valeurs moyennes de Cpa,Cpv et Lv
# en faisant le calcul de l'air moyen (i.e. au centre du diagramme
# le calcul s'arrête à 1e-6 de différence sur Cpa,Cpv et Lv
    crit=True
    while crit:
        Ain.car=Ma.array(Ain.car,copy=0,mask=maskin) 
        Aout.car=Ma.array(Aout.car,copy=0,mask=maskout) 
        Fill_Air(Ain)
        Fill_Air(Aout)
        Amoy=air()
        Amoy.definir(["H","na"],[mean([Ain.car[0],Aout.car[0]]),
                                    mean([Ain.car[1],Aout.car[1]])])
        Fill_Air(Amoy)
        Tmoy=Amoy.car[2]
        Cpa_new=Cpas(C2K(Tmoy))
        Cpv_new=Cpvap(C2K(Tmoy))
        Lv_new=ChLv(C2K(Tmoy))
        crit=((Cpa-Cpa_new)**2+(Cpv-Cpv_new)**2+(Lv-Lv_new)**2)>1e-6
        Cpa=Cpa_new
        Cpv=Cpv_new
        Lv=Lv_new
    if fignum>0 :
        Ratiofig=21.0/29.7
        fig_hauteur=6 #inches
        fig_longueur=fig_hauteur/Ratiofig
        figure(fignum,figsize=(fig_longueur,fig_hauteur))
        lineprops = dict(linewidth=0.5, color='gray', linestyle='-',antialiased=True)
        lineprops_HR = dict(linewidth=1, color='blue', linestyle='-',antialiased=True)
        LH,Lna,LT,LHR,Lair=Return_levels_HnaTHR(Ain,Aout)
        x1,y1=Transf_xy(Lna[0],Lair[0].car[0])
        x3,y3=Transf_xy(Lna[-1],Lair[2].car[0])
        RatioAxes=(y3-y1)/(x3-x1)
        return LH,Lna,LT,LHR,Lair

def Trace_diag(LH,Lna,LT,LHR,Lair,fignum) :
    figure(fignum)
    if isinteractive() :
        ioff()  
    lna_H,lH_H=Return_Hs(LH,Lna)
    Plotlines(lna_H,lH_H,lineprops)
    lna_T,lH_T=Return_Ts(LT,Lna)
    Plotlines(lna_T,lH_T,lineprops)
    lna_HR,lH_HR=Return_HRs(Lna,LHR)
    PlotHR(lna_HR,lH_HR,LHR,lineprops_HR)
    Ttick=Transf_xy(Lna[0],array(lH_T[0:-1:2]))[1]
    Tlabel=map(str,LT)
    yticks(Ttick,Tlabel)
    xlabel(r'$\rm{Teneur\ en\ eau\ (na)}\ \left[g_{H_2O}/g_{AS}\right]$')
    ylabel(r'$\rm{Temp\'erature\ (T)}\ \left[^{\circ}C\right]$')
    title("Diagramme de l'air humide")
    if not(isinteractive()) :
        ion()
    axis([x1,x3,y1,y3])
    setp(gca(),autoscale_on=False)
    draw()

#    ____________
#___/ Air humide \_____________________________________________________________
# Diverses fonctions permettant de calculer les paramètres thermophysiques d'un 
# air humide 

    
def rhoah(Aw):
    """Calcul de la masse volumique d'un air humide "Aw" à pression "p" (en Pa)
    le résultat est en kg/m^3 """
    na=Aw.car[1]
    T=Aw.car[2]
    Peau=na/(na+RatioM)*Ptot
    Pas=Ptot-Peau
    return 1/C2K(T)/Rgaz*(Peau*MH2O+Pas*Mair)/1000.0

# Titre liquide de l'air humide
def XI(Aw):
    """ Calcul du titre liquide (masse d'eau sur masse totale) d'un
    air humide Aw  """
    return Aw.car[1]/(1+Aw.car[1])

# Capacité thermique
def Cpah(Aw):
    """ Calcul de la capacité thermique massique d'un air humide"""
    # Titre liquide de l'air humide.
    T=C2K(Aw.car[2])
    x=XI(Aw)
    return Cpvap(T)*x+Cpas(T)*(1-x)
    
# Viscosité dynamique
def muah(Aw):
    """ Calcul de la viscosité dynamique de l'air humide """
    # Titre liquide de l'air humide.
    T=C2K(Aw.car[2])
    x=XI(Aw)
    return muvap(T)*x+muas(T)*(1-x)

# Viscosité cinématique
def nuah(Aw):
    return muah(Aw)/rhoah(Aw)

# Conductivité thermique
def lambdaah(Aw):
    """ Calcul de la conductivité thermique de l'air humide """
    # Titre liquide de l'air humide.
    T=C2K(Aw.car[2])
    x=XI(Aw)
    return lambdavap(T)*x+lambdaas(T)*(1-x)

# Nombre de Prandtl
def Prah(Aw):
    """ Calcul du nombre de Prandtl de l'air humide """
    # Titre liquide de l'air humide.
    T=C2K(Aw.car[2])
    x=XI(Aw)
    return Prvap(T)*x+Pras(T)*(1-x)

if __name__=='__main__' :
    A_bg=air()
    A_bg.definir(["T","na"],[5.0,0.001])
    A_hd=air()
    A_hd.definir(["T","na"],[110.0,0.2])
    LH,Lna,LT,LHR,Lair=Init_diag(A_bg,A_hd,1)
    Trace_diag(LH,Lna,LT,LHR,Lair,1)
    
