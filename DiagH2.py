# -*- coding: iso-8859-15 -*-
# Program permettant de tracer un diagramme de Mollier
# pour de l'hydrogène.
# 
#======================================================

from pylab import draw,figure,plot,isinteractive,ion,ioff,axis,xlabel,ylabel,xticks,yticks,text,title,show,setp,gca
from numpy import array,where,mean,floor,ceil,arange,unique,size,exp,arctan2,pi,select,linspace,logspace,log10
from scipy.optimize import brentq
from string import join
import MA as Ma
import copy
from proprietes_fluides import *
try: import psyco; psyco.full() 
except: pass

class h2 :
	"""Classe qui permet de définir un h2 il est caractérisé par 
	* h2.car : les 6 caractéristiques [H,na,T,HR,Tr,Th] sous forme d'un 
	tableau masqué T,Tr,Th en °C, na en g eau / g h2 sec, H en J/kg h2 sec 
	* h2.status : une string
		* indefini
		* defini
		* brouillard
		* saturation
	* h2.header : les entêtes i.e. list ["H","na","T","HR","Tr","Th"] 
	======================================
	Fonctions possibles
	* h2.definir permet de définir l'h2 à l'aide de 2 coordonnées 
	exemple h2.definir(["T","na"],[20,0.015])
	le statut passe alors à défini
	* h2.ou_non_defini permet de savoir les valeurs definies 
	exemple Tot,index=h2.ou_non_defini()
		* Tot le nombre total de car. non définies
		* index les index des valeurs non définies
	* h2.ou_defini idem ci-dessus
	exemple Tot,index=h2.ou_defini()
	* h2.affiche()
	exemple print h2.affiche() -> [20000 0.015 -- 0.6 -- 17.5]
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
		idx=where(self.car.mask()==[1,1,1,1,1,1])
		return [6.0-self.car.count(),idx]

	def ou_defini(self) :
		idx=where(self.car.mask()==[0,0,0,0,0,0])
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
	Func_calTh(T,H,HR) = Cph*T*Ptot+(H+Lv*RatioM)*HR*Pvap(T)+...
	(RatioM*Cpv-Cph)*T*HR*Pvap(T)-H*Ptot
	"""
	return Cph*T*Ptot+(H+Lv*RatioM)*HR*Pvap(T)+(RatioM*Cpv-Cph)*T*HR*Pvap(T)-H*Ptot

def Fill_from_Hna(h2):
	"""Permet de calculer les caractéristiques d'un h2 à partir des 
	coordonnées de référence H et na le statut de l'h2 passe alors à 
	brouillard ou saturation si besoin"""
	if h2.status=="indefini" :
		print "h2 indefini calcul impossible"
	else :
		nondef=h2.ou_non_defini() 
		H=h2.car[0]
		na=h2.car[1]
		T1,T2=Tminopt,h2.car[2]
		T3,T4=h2.car[4],h2.car[2]
		for i in nondef[1][0] :
			if i==2 :
				h2.car[i]=(H-na*Lv)/(Cph+Cpv*na)
				T2,T4=h2.car[i],h2.car[i]
			elif i==3 : 
				h2.car[i]=float(na/(na+RatioM)*Ptot/Pvap(h2.car[2]))
				if h2.car[i]>1.0 :
					h2.status="brouillard"
					#print "h2 dans une zone de brouillard"
					T1,T2=h2.car[2],Tmaxcrit
			elif i==4 :
				if h2.car[3]==1.0 :
					h2.car[i]=h2.car[2]
					h2.status="saturation"
				else :
					#print T1,T2,h2.affiche()
					h2.car[i]=brentq(Func_calTR,T1,T2,args=(h2.car[1],1.0))
					if h2.status=="brouillard" : 
						T3,T4=h2.car[2],h2.car[i]
					else :
						T3,T4=h2.car[i],h2.car[2]
			elif i==5 :
				if h2.car[3]==1.0 :
					h2.car[i]=h2.car[2]
					h2.status="saturation"
				else :
					print T3,T4,h2.car
					h2.car[i]=brentq(Func_calTh,T3,T4,args=(h2.car[0],1.0))
		if h2.status=="brouillard" :
			pass
		else :
			h2.status="defini"
	return h2

def Fill_h2(h2):
	"""Permet de calculer les caractéristiques d'un h2 à partir de 2 coordonnées 
	par le jeu des relations on se ramème au cas H et na"""
	estdef=h2.ou_defini()
	#print "La Pression est de :",Ptot
	while (estdef[0]<6) :
		Prop_input=estdef[1][0][0:2]
		# on verifie si on est dans une zone possible
		# verif sur H
		if (h2.car.mask()[0]==0) :
			if (h2.car[0]<Hmincrit) :
				h2.car[0]=Hmincrit
				Fill_h2(h2)
				H=h2.car[0]
			else :
				H=h2.car[0]
		# verif sur na
		if (h2.car.mask()[1]==0) :
			if h2.car[1]<namincrit :
				h2.car[1]=namincrit
				Fill_h2(h2)
				na=h2.car[1]
			else :
				na=h2.car[1]
		# verif sur T
		if (h2.car.mask()[2]==0) :
			if h2.car[2]<Tmincrit :
				h2.car[2]=Tmincrit
				Fill_h2(h2)
				T=h2.car[2]
			else :
				T=h2.car[2]
		# pas de verif sur HR car > 1 est possible : brouillard
		if (h2.car.mask()[3]==0) :	
			HR=h2.car[3]
		# verif sur TR
		if (h2.car.mask()[4]==0) :	
			if h2.car[4]<Tmincrit :
				h2.car[4]=Tmincrit
				Fill_h2(h2)
				TR=h2.car[4]
			else :
				TR=h2.car[4]
		#verif sur Th
		if (h2.car.mask()[5]==0) :	
			if h2.car[5]<Tmincrit :
				h2.car[5]=Tmincrit
				Fill_h2(h2)
				Th=h2.car[5]
			else :
				Th=h2.car[5]
		
		# From H and na
		if (Prop_input==[0,1]).all() :
			Fill_from_Hna(h2)
		# From H and T
		elif (Prop_input==[0,2]).all() :
			h2.car[1]=float((H-Cph*T)/(Lv + Cpv*T))
			Fill_from_Hna(h2)
		# From H and HR
		elif (Prop_input==[0,3]).all() :
			h2.car[2]=brentq(Func_calTh,Tminopt,Tmaxcrit,args=(H,HR))
		# From H and TR
		elif (Prop_input==[0,4]).all() : 
			h2.car[1]=float(RatioM*Pvap(TR)/(Ptot-Pvap(TR)))
		# From H and Th
		elif (Prop_input==[0,5]).all() : 
			print "Operation impossible h2 indefini car H et Th"
			break
		# From na and T
		elif (Prop_input==[1,2]).all() : 
			h2.car[0]=float(Cph*T+na*(Lv+Cpv*T))
		# From na and HR
		elif (Prop_input==[1,3]).all() : 
			h2.car[2]=brentq(Func_calTR,Tminopt,Tmaxcrit,args=(na,HR))
		# From na and TR
		elif (Prop_input==[1,4]).all() : 
			print "Operation impossible h2 indefini car na et TR"
			break
		# From na and Th
		elif (Prop_input==[1,5]).all() : 
			h2.car[0] = float(Cph*Th +RatioM*Pvap(Th)/(Ptot-Pvap(Th))*(Lv + Cpv*Th))
		# From T and HR
		elif (Prop_input==[2,3]).all() : 
			h2.car[1]=float(RatioM*Pvap(T)*HR/(Ptot-Pvap(T)*HR))
		# From T and TR
		elif (Prop_input==[2,4]).all() : 
			h2.car[1]=float(RatioM*Pvap(TR)/(Ptot-Pvap(TR)))
		# From T and Th
		elif (Prop_input==[2,5]).all() : 
			h2.car[0] =float( Cph*Th +RatioM*Pvap(Th)/(Ptot-Pvap(Th))*(Lv + Cpv*Th))
		# From HR and TR
		elif (Prop_input==[3,4]).all() : 
			h2.car[1]=float(RatioM*Pvap(TR)/(Ptot-Pvap(TR)))
		# From HR and Th
		elif (Prop_input==[3,5]).all() : 
			h2.car[0] = float(Cph*Th +RatioM*Pvap(Th)/(Ptot-Pvap(Th))*(Lv + Cpv*Th))
		# From TR and Th
		elif (Prop_input==[4,5]).all() :
			h2.car[0] = float(Cph*Th +RatioM*Pvap(Th)/(Ptot-Pvap(Th))*(Lv + Cpv*Th))
		else : pass
		estdef=h2.ou_defini()
	return h2

def Return_BBox(h2_bg,h2_hd) :
	"""Fonction qui retourne les caractéristiques des h2s aux 
	4 coins du diagramme à partir de l'h2 en bas à gauche et 
	celui en haut à droite"""
	h2_hg=h2()
	h2_hg.definir(["na","H"],[h2_bg.car[1],h2_hd.car[0]-(h2_hd.car[1]-h2_bg.car[1])*Lv])
	Fill_h2(h2_hg)
	h2_bd=h2()	
	h2_bd.definir(["na","H"],[h2_hd.car[1],h2_bg.car[0]-(h2_bg.car[1]-h2_hd.car[1])*Lv])
	Fill_h2(h2_bd)
	h21=copy.deepcopy(h2_bg)
	h23=copy.deepcopy(h2_hd)
	return h21,h2_bd,h23,h2_hg

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

def Return_levels_HnaTHR(h2in,h2out) :
	"""Fonction qui renvoie les niveaux de H,na,T et HR à afficher"""
	h21,h22,h23,h24=Return_BBox(h2in,h2out)  
	namin=max(floor(h21.car[1]*1000)/1000,namincrit)
	namax=ceil(h22.car[1]*1000)/1000
	#print namin,namax
	#print "h21",h21.affiche()
	#print "h22",h22.affiche()
	#print "h23",h23.affiche()
	#print "h24",h24.affiche()
	
	if floor(h21.car[2])<floor(h22.car[2]) :
		Tmin=floor(h21.car[2])
		h21.definir(["na","T"],[namin,Tmin])
		Fill_h2(h21)
		h22.definir(["na","H"],[namax,h21.car[0]-(namin-namax)*Lv])
		Fill_h2(h22)
		#print "Cas 1 T1<T2 Tmin=",Tmin		
	else :
		Tmin=floor(h22.car[2])
		h22.definir(["na","T"],[namax,Tmin])
		Fill_h2(h22)
		h21.definir(["na","H"],[namin,h22.car[0]-(namax-namin)*Lv])
		Fill_h2(h21)	
		#print "Cas 2 T1>T2 Tmin=",Tmin
		
	if ceil(h23.car[2])<ceil(h24.car[2]) :
		Tmax=ceil(h24.car[2])
		h24.definir(["na","T"],[namin,Tmax])
		Fill_h2(h24)
		h23.definir(["na","H"],[namax,h24.car[0]-(namin-namax)*Lv])
		Fill_h2(h23)
	else :
		Tmax=ceil(h23.car[2])
		h23.definir(["na","T"],[namax,Tmax])
		Fill_h2(h23)
		h24.definir(["na","H"],[namax,h23.car[0]-(namax-namin)*Lv])
		Fill_h2(h24)
	
	HRmin=floor(h24.car[3]*100)/100
	HRmax=ceil(min(1.0,h22.car[3])*100)/100
	Hmin=floor(h21.car[0]*1e3)/1e3
	Hmax=ceil(h23.car[0]*1e3)/1e3	
	na_base=arange(0.0,2.0,0.0005)
	HR_base=arange(0.0,1.001,0.01)
	T_base=arange(Tmincrit,Tmaxcrit+0.01,0.1)
	H_base=arange(-411000,5e6+1.0,1e3)
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
	#print "h21",h21.affiche()
	#print "h22",h22.affiche()
	#print "h23",h23.affiche()
	#print "h24",h24.affiche()
	return list_H,list_na,list_T,list_HR,[h21,h22,h23,h24]

def Transf_xy(x,y) :
	"""Effectue une tranformation d'axe
	de na,H en x,y"""
	y=y-x*Lv
	return x,y

def Return_Hs(lh,lna) :
	"""Permet de calculer des isoH avec
	* h2_1 et h2_2 en entree
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
	* h2_1 et h2_2 en entree
	* n = nombre d'isoT"""
	list_H=[]
	list_na=[]
	h2tmp=h2()
	for i in lt :
		h2tmp.definir(["T","na"],[i,lna[0]])
		Fill_h2(h2tmp)
		list_H.append(h2tmp.car[0])
		h2tmp.definir(["T","na"],[i,lna[-1]])
		Fill_h2(h2tmp)
		list_H.append(h2tmp.car[0])	
	list_na.append(lna[0])
	list_na.append(lna[-1])
	return list_na,list_H
	
def Return_HRs(lna,lhr) :
	"""Permet de calculer des isoHR avec
	* """
	list_H=[]
	metalist_H=[]
	list_na=[]
	h2tmp=h2()
	namin=lna[0]
	namax=lna[-1]
	A=10.
	e=exp(A)
	C=(namax-namin)/(e-1)
	Vecna=C*(exp(A*(lna-namin)/(namax-namin))-1)+namin
	#Vecna=logspace(log10(namin),log10(namax),num=100,endpoint=True)
	for i in lhr :
		for j in Vecna :
			h2tmp.definir(["HR","na"],[float(i),float(j)])
			Fill_h2(h2tmp)
			#print h2tmp.affiche()
			list_H.append(h2tmp.car[0])
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

def Init_diag(h2in,h2out,fignum):
	""" Procédure d'initialisation du diagramme et qui permet de fixer 
	les limites du diagramme. Par exemple, on fixe h2insi dans une variable 
	globale la pression totale. Voilà un exemple d'utilisation :
	
	n=1
	h2_bg=Diag.h2()
	h2_bg.definir(["T","na"],[20.0,0.00001])
	h2_hd=Diag.h2()
	h2_hd.definir(["T","na"],[80.0,0.6])
	LH,Lna,LT,LHR,Lh2=Diag.Init_diag(h2_bg,h2_hd,n)
	
	avec h2_bg l'h2 en bas à gauche et h2_hd l'h2 en haut à droite
	le dernier paramètre ici n=1 spécifie le numéro de la figure matplotlib. 
	LH,Lna,LT,LHR,Lh2 sont les sorties et correpondent aux niveaux sur
	les différentes coordonnées. Lh2 et une liste de taille 4 qui contient 
	les h2s des 4 coins du diagramme. 
	
	si l'on met n=-1 alors on ne fait que des calculs sur des h2s humides
	et l'appel se fait h2insi :
	
	Diag.Init_diag(h2_bg,h2_hd,-1)
	"""
	global Ptot,RatioM,Tmincrit,Cph,Cpv,Lv,namincrit,Hmincrit,Tmaxcrit,Tminopt,Rgaz,MH2O,MH2
	global Ratiofig,fig_hauteur,fig_longueur,RatioAxes,lineprops,lineprops_HR,x1,y1,x3,y3
	Rgaz=8.314472
	Ptot=101325.0
	MO=15.9994
	MH=1.00794
	MH2=2.*MH
	MH2O=2*MH+MO
	RatioM=MH2O/MH2
	Tmincrit=-23.0
	Cph=Cph2(C2K(Tmincrit))
	Cpv=Cpvap(C2K(0.1))
	Lv=ChLv(C2K(0.1))
	namincrit=RatioM*Pvap(Tmincrit)/(Ptot-Pvap(Tmincrit))
	Hmincrit=Cph*Tmincrit+namincrit*(Lv+Cpv*Tmincrit)
	Tmaxcrit=125.0
	Tminopt=-30.0
	Tinit=20.0
	Cph=Cph2(C2K(Tinit))
	Cpv=Cpvap(C2K(Tinit))
	Lv=ChLv(C2K(Tinit))
	maskin=Ma.getmask(h2in.car)
	maskout=Ma.getmask(h2out.car)
	
# Il s'agit ici de fixer les valeurs moyennes de Cph,Cpv et Lv
# en faisant le calcul de l'h2 moyen (i.e. au centre du diagramme
# le calcul s'arrête à 1e-6 de différence sur Cph,Cpv et Lv
	crit=True
	while crit:
		h2in.car=Ma.array(h2in.car,copy=0,mask=maskin) 
		h2out.car=Ma.array(h2out.car,copy=0,mask=maskout) 
		Fill_h2(h2in)
		Fill_h2(h2out)
		h2moy=h2()
		h2moy.definir(["H","na"],[mean([h2in.car[0],h2out.car[0]]),
						            mean([h2in.car[1],h2out.car[1]])])
		Fill_h2(h2moy)
		Tmoy=h2moy.car[2]
		Cph_new=Cph2(C2K(Tmoy))
		Cpv_new=Cpvap(C2K(Tmoy))
		Lv_new=ChLv(C2K(Tmoy))
		crit=((Cph-Cph_new)**2+(Cpv-Cpv_new)**2+(Lv-Lv_new)**2)>1e-6
		Cph=Cph_new
		Cpv=Cpv_new
		Lv=Lv_new
	if fignum>0 :
		Ratiofig=21.0/29.7
		fig_hauteur=6 #inches
		fig_longueur=fig_hauteur/Ratiofig
		figure(fignum,figsize=(fig_longueur,fig_hauteur))
		lineprops = dict(linewidth=0.5, color='gray', linestyle='-',antialiased=True)
		lineprops_HR = dict(linewidth=1, color='blue', linestyle='-',antialiased=True)
		LH,Lna,LT,LHR,Lh2=Return_levels_HnaTHR(h2in,h2out)
		x1,y1=Transf_xy(Lna[0],Lh2[0].car[0])
		x3,y3=Transf_xy(Lna[-1],Lh2[2].car[0])
		RatioAxes=(y3-y1)/(x3-x1)
		return LH,Lna,LT,LHR,Lh2

def Trace_diag(LH,Lna,LT,LHR,Lh2,fignum) :
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
	xlabel(r'$\rm{Teneur\ en\ eau\ (na)}\ \left[g_{H_2O}/g_{H_2}\right]$')
	ylabel(r'$\rm{Temp\'erature\ (T)}\ \left[^{\circ}C\right]$')
	title(u"Diagramme de l'hydrogène humide")
	if not(isinteractive()) :
		ion()
	axis([x1,x3,y1,y3])
	setp(gca(),autoscale_on=False)
	draw()

#    ___________
#___/ h2 humide \_____________________________________________________________
# Diverses fonctions permettant de calculer les paramètres thermophysiques d'un 
# h2 humide 

	
def rhoh2h(h2w):
	"""Calcul de la masse volumique d'un h2 humide "hw" à pression "p" (en Pa)
	le résultat est en kg/m^3 """
	na=h2w.car[1]
	T=h2w.car[2]
	Peau=na/(na+RatioM)*Ptot
	Pas=Ptot-Peau
	return 1./C2K(T)/Rgaz*(Peau*MH2O+Pas*MH2)/1000.0

# Titre liquide de l'h2 humide
def XIh2(h2w):
	""" Calcul du titre liquide (masse d'eau sur masse totale) d'un
	h2 humide Aw  """
	return h2w.car[1]/(1.+h2w.car[1])

# Capacité thermique
def Cph2h(h2w):
	""" Calcul de la capacité thermique massique d'un h2 humide"""
	# Titre liquide de l'air humide.
	T=C2K(h2w.car[2])
	x=XI(h2w)
	return Cpvap(T)*x+Cph2(T)*(1.-x)
	
# Viscosité dynamique
def muh2h(h2w):
	""" Calcul de la viscosité dynamique de l'h2 humide """
	# Titre liquide de l'air humide.
	T=C2K(h2w.car[2])
	x=XI(h2w)
	return muvap(T)*x+muh2(T)*(1.-x)

# Viscosité cinématique
def nuh2h(h2w):
	return muh2h(h2w)/rhoh2h(h2w)

# Conductivité thermique
def lambdah2h(h2w):
	""" Calcul de la conductivité thermique de l'h2 humide """
	# Titre liquide de l'air humide.
	T=C2K(h2w.car[2])
	x=XIh2(h2w)
	return lambdavap(T)*x+lambdah2(T)*(1.-x)

# Nombre de Prandtl
def Prh2h(h2w):
	""" Calcul du nombre de Prandtl de l'h2 humide """
	# Titre liquide de l'air humide.
	T=C2K(h2w.car[2])
	x=XIh2(h2w)
	return Prvap(T)*x+Prh2(T)*(1.-x)


