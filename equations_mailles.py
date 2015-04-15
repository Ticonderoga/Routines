# -*- coding: iso-8859-15 -*-

# -------------------------------------------------------------------
# |             Sous programmes de resolution d'equations de mailles            |
# -------------------------------------------------------------------

#    ________________________
#___/ Matrices tridiagonales \________________________________________

def eqmailles_balay_droite(am,a,ap,b):
        # Resolution de l'equation de mailles representee par le systeme
        # lineaire tridiagonale avec "am", "a" et "ap" respectivement les
        # diagonales inferieure, principale et superieure de la matrice.
        # "b" est le second membre du systeme. 
        # ([1], Chapitre 2, §1, pp 77)
        # Les trois diagonales "am" (de taille N-1),"ap" (de taille N-1)  
        # et "a" (de taille N) sont entrees dans le sous programme sous 
        # formes de listes afin de minimiser l'espace memoire. Le second 
        # membre de l'equation est appele "b" (il est de dimensions N,M) 
        # et la solution du systeme "x" (de dimension N,M). M étant le
        # nombre de seconds membres que l'on traite avec une même matrice.
        # On verifie tout d'abord que les dimensions des differentes
        # diagonales sont compatibles entre elles.
        if ( ( len(am) != len(ap) ) or ( len(am) != (len(a)-1) ) ):
                # Message d'erreur et arrêt dans le cas contraire.
                print "eqmailles_balay_droite : (Erreur) Probleme de taille"
                print "dans les arguments d'entree. Certaines des diagonales"
                print "ont des tailles incompatibles entre elles."
        else:
                # On teste ensuite la dimension du second membre.
                if ( (type(b[0]) == float) or (type(b[0]) == int) ):
                        # Dans le cas où le second membre est un vecteur simple,
                        # on appelle le sous programme concerne.
                        return eqmailles_balay_droite_scal(am,a,ap,b)
                else:
                        # Si par contre il s'agit d'un ensemble de vecteurs, le sous
                        # programme est different.
                        return eqmailles_balay_droite_vect(am,a,ap,b)
#---------------------------------------------------------------------
def eqmailles_balay_droite_scal(am,a,ap,b):
    # Application de la methode de balayage à droite pour les systemes 
    # tridiagonaux dont le second membre est un vecteur simple.
    # Cette fonction est utilisee par :
    # - eqmailles_balay_droite
    # On effectue un test de compatibilité sur les dimensions de de la
    # matrice et du second membre.
    if ( len(a) != len(b) ):
        # Message d'erreur et arrêt :
        print "eqmailles_balay_droite_scal : (Erreur) Probleme de taille"
        print "dans les arguments d'entree. La taille de la matrice et du"
        print "second membre ne sont pas compatibles entre elles."
    else:
        # Dimension du probleme :
        Nx=len(a)
        # On commence par calculer les premiers coefficients de balayage
        alpha=coeff_balayage_01_dr(am,a,ap)
        # Puis les seconds, qui dependent notamment du second membre
        beta=coeff_balayage_02_dr(am,a,ap,b,alpha)
        # Et on calcule enfin les termes du vecteur recherche.
        return solution_balayage_dr(alpha,beta)

#---------------------------------------------------------------------
def eqmailles_balay_droite_vect(am,a,ap,b):
        # Application de la methode de balayage à droite pour les systemes 
        # tridiagonaux dont le second membre est une suite de vecteurs.

        # Cette fonction est utilisee par :
        # - eqmailles_balay_droite

        # Dimension du probleme :
        Nx=len(a)
        # On commence par calculer les premiers coefficients de balayage
        alpha=coeff_balayage_01_dr(am,a,ap)
        # On determine ensuite les dimensions du second membre.
        d1=len(b)
        d2=len(b[0])
        # On effectue ensuite une serie de test sur les dimensions de b
        if ((d1!=len(a)) and (d2!=len(a))):
                # Si aucune des dimensions ne convient à celle du systeme, on
                # envoie un message d'erreur.
                print "eqmailles_balay_droite : (Erreur) Probleme de"
                print "taille dans les arguments d'entree. Le second"
                print "membre du systeme est d'une taille incompatible"
                print "avec la matrice."
        elif ( (d1==d2) and (d1==len(a)) ):
                # Si au contraire les deux conviennent, il y a une
                # ambiguite sur la maniere dont sont assembles les
                # vecteurs et un avertissement apparâit.
                print "eqmailles_balay_droite : (Avertissement)"
                print "Possible ambiguite dans la taille du second"
                print "membre du systeme, verifiez soigneusement vos"
                print "arguments d'entree."
        # Si les seconds membres sont presentes de "maniere
        # transposee", on modifie leurs dispositions.
        if ((d1==len(a)) and (d2!=len(a))):
                bb=b
                for i in range(0,d1):
                        for j in range(0,d2):
                                b[i][j]=bb[j][i]
        # On effectue alors une serie de resolutions en n'appelant
        # que la seconde serie de coefficients de balayages, pour
        # chaque vecteur du second membre.
        # On crée d'abord la matrice des solutions
        sol=[[float(0)]*Nx]*d1
        # Puis on lance la résolution pour chacun des seconds membres
        for k in range(0,d1):
                # On appelle les coefficients de balayage correspondants
                # à chaque second membre.
                beta=coeff_balayage_02_dr(am,a,ap,b[k],alpha)
                # Initialisation d'une liste temporaire.
                sol2=[float(0)]*Nx
                # En on commence le calcul des solutions par le dernier 
                # terme de chaque vecteur.
                sol2[Nx-1]=beta[Nx-1]
                # Puis on traite les termes précédents
                for i in range(1,Nx):
                        sol2[Nx-i-1]=beta[Nx-i-1]+alpha[Nx-i-1]*sol2[Nx-i]
                sol[k]=sol2
        return sol
#---------------------------------------------------------------------
def coeff_balayage_01_dr(am,a,ap):
        # Programme de calcul du premier jeux de coefficients de balayage,
        # ceux qui ne dependent que de la matrice tridiagonale. La methode
        # de balayage consideree est celle "à droite" et les trois 
        # arguments d'entree sont "am", "a" et "ap" respectivement les
        # diagonales inferieure, principale et superieure.
        # La sortie est la liste de ces coefficients de balayage, de taille
        # N-1 si la matrice est de taille N.
        if ( len(am) != len(ap) ) or ( len(am) != (len(a)-1) ):
                print "coeff_balayage_01_dr : (Erreur) Probleme de taille dans"
                print "les arguments d'entree"
        else:
                # On calcule ensuite ces coefficients de balayage, en 
                # commençant par les initialiser à 0.
                alpha=[float(0)]*(len(a)-1)
                # Puis on calcule le premier coefficient.
                alpha[0]=-float(ap[0]/a[0])
                # Et les suivants
                for i in range(1,len(a)-1):
                        alpha[i]=-float(ap[i]/(a[i]+(am[i-1]*alpha[i-1])))
                return alpha
#
def coeff_balayage_02_dr(am,a,ap,b,alpha):
        # Calcul de la deuxieme serie de coefficients de balayage, ceux qui
        # dependent aussi du second membre et des premiers coefficients.
        # L'execution de ce programme ne pouvant se faire sans celui qui
        # calcule les premiers coefficients, on se passera ici d'appliquer
        # les mêmes procedures de tests.
        beta=[float(0)]*len(a)
        beta[0]=b[0]/a[0]
        for i in range(1,len(a)):
                beta[i]=(b[i]-am[i-1]*beta[i-1])/(a[i]+am[i-1]*alpha[i-1])
        return beta
#
def solution_balayage_dr(alpha,beta):
    # Calcul de la solution d'un système linéaire tridiagonal à partir
    # seulement des coefficients de balayage à droite.
    
    # Cette fonction est utilisee par :
    # - eqmailles_balay_droite_scal
    
    # Dimension du système.
    Nx=len(beta)
    # Initialisation de la solution
    sol=[float(0)]*Nx
    # Et calcul de cette solution à partir des coefficients de balayage
    sol[Nx-1]=beta[Nx-1]
    for i in range(1,Nx):
        sol[Nx-i-1]=beta[Nx-i-1]+alpha[Nx-i-1]*sol[Nx-i]
    return sol

# References :
# ------------

# [1] Samarski, A. and E. Nikolaïev, Methodes de Resolution des  
#     Equations de Mailles, ed. Mir. 1981, Moscou.
