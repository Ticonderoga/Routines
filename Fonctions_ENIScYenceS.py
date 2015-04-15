# -*- coding: iso-8859-15 -*-

import numpy as Nu
import scipy.sparse as Sp
from GVar_Stirling import *


def ParseResults(P,R,Bool) :
	return [Sp.csr_matrix(Nu.where(Bool,P,0.)),Sp.csr_matrix(Nu.where(Bool,R,0.))]
	

def Pow_Eff(teta,DT,epsi,gama,alpha,Kh,Kc,Kl,Kreg,R,Tc,Th) :
	""" Calcul de la puissance et du rendement d'un moteur Stirling
	C'est une relation vectorielle et les variables :
		* teta et DT sont des vecteurs.
		* les autres parmètres sont des scalaires
	Exemple d'utillisation :
	[mat_P,mat_R]=Pow_Eff(vec_teta,vec_DT,epsi,gama,alpha,Kh,Kc,Kl,Kreg,R,Tc,Th)
		* mat_P et mat_R sont des matrices creuses """
	teta=Nu.atleast_1d(Nu.asarray(teta,dtype=float))
	DT=Nu.atleast_1d(Nu.asarray(DT,dtype=float))
	M_T=Nu.tile(teta,(DT.size,1))
	M_DT=Nu.transpose(Nu.tile(DT,(teta.size,1)))
	C1=(gama-1.)*(M_T-R)/(M_T-1.)*Nu.log(epsi)
	C2=Nu.log((M_T-alpha*(M_T-1.))*((1.+alpha*(M_T-1.))/M_T)**R)/(M_T-1.)
	C3=(Nu.log(epsi)*(gama-1.)+(1.-alpha)*(M_T-1))/(Kc*(M_DT-Tc*(M_T-1.)))
	C4=(M_T*Nu.log(epsi)*(gama-1.)+(1.-alpha)*(M_T-1))/(Kh*(Th*(M_T-1.)-M_T*M_DT))
	C5=2./Kreg/M_DT
	A1=C1+C2
	A2=C3+C4+C5
	Pow=Nu.squeeze(A1/A2/1000.) # Puissance en kW
	Eff=Nu.squeeze(A1/(M_T/(M_T-1.)*(gama-1.)*Nu.log(epsi)+(1.-alpha)+Kl*(Th-Tc)*A2))
	if Nu.size(Eff)==1 and Nu.size(Pow)==1:
		return [Pow.astype('float'),Eff.astype('float')]
	else :
		X,Y=Nu.meshgrid(teta,DT)
		Bool1=(Y<Th*(X-1.)/X)
		Bool2=Y>Tc*(X-1.)
		Bool=Bool1*Bool2
		return ParseResults(Pow,Eff,Bool)

