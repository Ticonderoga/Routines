# -*- coding: iso-8859-15 -*-

# -------------------------------------------------------------------
# |    Sous programmes divers utilisés dans différents types de      |
# |    calculs.                                                      |
# -------------------------------------------------------------------

from os import system
import numpy as np
import pylab as Py
try: import psyco; psyco.full() 
except: pass

def maillage(x,style,parameter,sens):
    """ Fonction permettant de remailler à la paroi en 1D
        soit : 
        x / le vecteur de valeurs initiales 
        style / la fonction de maillage (puissance ou exponentielle)
                style='power' ou style='exponential'
        parameter / le paramètre permettant de jouer sur la courbure
                   attention on a parameter > 1
        sens / permet de savoir si la paroi est à droite ou à gauche
                sens='right' ou sens='left'
        exemple avec un profil parabolique:
            x=np.linspace(45,70,6)
            # i.e. x= array([ 45.,  50.,  55.,  60.,  65.,  70.])
            y=maillage(x,'power',2,'right')
            # i.e. y=array([ 45.,  54.,  61.,  66.,  69.,  70.])
            y=maillage(x,'power',2,'left')
            # i.e. y=array([ 45.,  46.,  49.,  54.,  61.,  70.])"""
    xmin=x.min()
    xmax=x.max()
    xs=(x-xmin)/(xmax-xmin)
    if parameter>=1. :
        if style=='power' :
            f=lambda x : x**parameter
        elif style=='exponential' :
            f=lambda x : (1.-np.exp(parameter*x))/(1.-np.exp(parameter))
        if sens=='left' :
            return np.sort(f(xs)*(xmax-xmin)+xmin)
        elif sens=='right' :
            return np.sort((1.-f(xs))*(xmax-xmin)+xmin)
    else :
        print("You should have parameter > 1")

def test_longueur_liste(liste,dim):
    # Fonction de test de longueur d'une liste de valeurs. "liste" désigne la
    # liste à tester et "dim" le nombre désignant la dimension. Le programme
    # renvoie la valeur 0 si le test et positif et 1 dans le cas contraire.

    # Bien que la tâche traitée soit élémentaire, une telle fonction est
    # parfois appelée dans un programme python lorsque le type d'un argument
    # (liste ou scalaire) n'est pas définit de manière univoque.

    # On vérifie que la liste à tester est bien une liste.
    if (type(liste) != list):
        print("test_longueur_liste (Erreur) : Le premier argument envoye dans")
        print("la fonction n'est pas une liste.")
        return 1
    # Que le nombre représentant la dimension est bien un nombre.
    elif (type(dim) == list):
        print("test_longueur_liste (Erreur) : Le deuxième argument envoye dans")
        print("la fonction n'est pas une valeur scalaire.")
        return 1
    else:
        # Que le nombre représentant la dimension est bien un nombre entier.
        if (type(dim) == float):
            print("test_longueur_liste (Avertissement) : La dimension de la")
            print("liste a tester est entree sous la forme d'un nombre")
            print("flottant -> conversion automatique en nombre entier.")
        # On effectue enfin le test en question.
        if (len(liste) == dim):
            return 0
        else:
            return 1

def comp_longueur_liste(liste1,liste2):
    # Comparaison des longueurs de deux listes différentes. Cette fonction
    # renvoie "1" si les deux listes sont de même taille et "0" sinon.
    if (len(liste1) != len(liste2)):
    # Message d'erreur si listes de tailles différentes.
        print("comp_longueur_liste : (Erreur) Les deux listes à manipuler ne")
        print("sont pas de la même dimension.")
        return 0
    else:
        return 1

def interp_polyn_Lagrange(xi,yi,x):
    """ Sous programme d'approximation de la valeur y = p(x) où "x" est connp
    et où p est le polynôme de Lagrange d'ordre N-1 passant par les N points 
    de coordonnées (xi,yi)."""
    # On teste tout d'abord la compatibilité des dimensions des listes "xi" et
    # "yi".
    if (comp_longueur_liste(xi,yi)==0):
    # Message d'erreur renvoyé par la fonction de test
        return
    else:
    # On lance ensuite le calcul :
        yn=0
        for i in range(len(xi)):
            # On parcours les abcisses xi en calculant à chaque fois le
            # npmérateur et le dénominateur du polynome".
            NL=1
            DL=1
            for j in range(len(xi)):
                if (i != j):
                    NL=NL*(x-xi[j])
                    DL=DL*(xi[i]-xi[j])
                # Et on renvoie le résultat.
            yn=yn+(NL/DL)*yi[i]
    return yn

#def bin(n):
#    """Convertit un nombre en binaire"""
#    res = ''
#    while n != 0: n, res = n >> 1, `n & 1` + res
#    return res
    

def test_io(StrFile,ncol):
    """Fonction permettant d'ouvrir un fichier ASCII comportant un
    certain nombre de colonnes :
    * StrFile : String = le nom du fichier i.e. path+file
    * ncol : Integer = le nombre de colonnes """
    inFile = file(StrFile, 'r')
    Tab= array([float(x) for x in inFile.read().split()])
    nrow=Tab.size/ncol
    Tab.shape = (nrow,ncol)
    return Tab

def multiplot(Tab,*args) :
    """Fonction permettant de tracer un tableau selon les colonnes et
    la 1ère colonne correspond à l'axe des abscisses :
    * Tab : Array = Tableau à tracer 
            1ère col = axe des abscisses
    * *args : Tuple = Les labels pour la légende
    Exemple :
    label_data=("label_1","label_2") 
    multiplot(Donnees_3_colonnes,*label_data)"""
    
    hold(True)
    if len(args)==0 :   
        for i in range(1,size(Tab,1)) :
            plot(Tab[:,0],Tab[:,i])
    else :
        if size(args)!=(size(Tab,1)-1) :
            print("Discordance entre les étiquettes et les colonnes")
            print("Labels : ",size(args))
            print("Colonnes : ", size(Tab,1)-1 )
        else :
            for i in range(1,size(Tab,1)) :
                plot(Tab[:,0],Tab[:,i],label=args[i-1])
            legend(loc="best")  
            leg = gca().get_legend()
            # all the text.Text instance in the legend
            ltext  = leg.get_texts()
            # all the lines.Line2D instance in the legend  
            llines = leg.get_lines() 
            # the legend text fontsize
            setp(ltext, fontsize='small')
            # the legend linewidth
            #setp(llines, linewidth=1.5)      
            # don't draw the legend frame
            leg.draw_frame(False)           
    

def plot_Elec(Tab,nplot,ncurve) :
    """ Fonction de tracé des données électriques d'une PAC. Cela 
    permet de choisir le nombre de subplot et le nombre 
    de courbes à tracer :
    * Tab : Array = Données électriques à tracer.
        Il est formé ainsi :
            1 col   = temps 
            2..n    = les tensions Ui (V)
            n+1:n+2 = U (V) et I (A) globaux
    * nplot : Tuple = (m,n) nombre de subplots ordonnés selon une 
        matrice (m,n)
    * ncurve : Integer = Le nombre de courbes à tracer.
    NB : Le dernier subplot comporte toujours U et I."""
    ncell=size(Tab,1)-3
    if (nplot[0]*nplot[1]-1)*ncurve>=ncell :
        time=Tab[:,0:1]
        nent=ncell//ncurve
        nrem=ncell%ncurve
        # on trace par paquet de ncurve
        for i in range(nent) :
            subplot(nplot[0],nplot[1],i+1)
            Tab2pl=hstack((time,Tab[:,1+ncurve*i:ncurve*(i+1)+1]))
            multiplot(Tab2pl)
        # on trace les dernieres courbes
        # puis on trace UI sur le dernier subplot
        subplot(nplot[0],nplot[1],nent+1)
        Tab2pl=hstack((time,Tab[:,nent*ncurve:nent*ncurve+nrem]))
        multiplot(Tab2pl)
        subplot(nplot[0],nplot[1],nent+2)
        Tab2pl=hstack((time,Tab[:,-2:]))
        multiplot(Tab2pl)
    else :
        print("pas suffisament de place pour tracer")
        print(locals())

def plot_Data(Tab,nplot,ncurve,labx,laby,*args) :
    """ Fonction de tracé des données d'une PAC. Cela 
    permet de choisir le nombre de subplot, le nombre 
    de courbes à tracer et les étiquettes des axes :
    * Tab : Array = Données à tracer.
        Il est formé ainsi :
            1 col   = temps (ms)
            2..n    = les données
    * nplot : Tuple = (m,n) nombre de subplots ordonnés selon une 
        matrice (m,n)
    * ncurve : Integer = Le nombre de courbes à tracer.
    labx : String = étiquette de l'axe des x
    laby : String = étiquette de l'axe des y
    * *args : Tuple = Les labels pour la légende Cf. multiplot"""
    ndata=size(Tab,1)-1
    if nplot[0]*nplot[1]*ncurve>=ndata :
        nent=ndata//ncurve
        nrem=ndata%ncurve
        time=Tab[:,0:1]
        for i in range(nent) :
            ax=subplot(nplot[0],nplot[1],i+1)
            Tab2pl=hstack((time,Tab[:,1+ncurve*i:ncurve*(i+1)+1]))
            argbis=args[ncurve*i:ncurve*(i+1)]
            multiplot(Tab2pl,*argbis)
            if ax.is_last_row() : xlabel(labx)
            if ax.is_first_col() : ylabel(laby)
        if nrem!=0 :
            subplot(nplot[0],nplot[1],nplot[0]*nplot[1])
            Tab2pl=hstack((time,Tab[:,-nrem:]))
            argbis=args[-nrem:]
            multiplot(Tab2pl,*argbis)
    else :
        print("pas suffisament de place pour tracer")
        print(locals())



def moypaq(Tab,n) :
    """ fait des moyennes selon n lignes
    Exemple : 
    A=array([[1, 2, 3, 4, 4, 6],
             [1, 2, 3, 1, 1, 1]])
    moypaq(A,2)=array([[1.,2.,3.,2.5,2.5,3.5]])
    """
    nlignes,reste=divmod(Tab.shape[0],n)
    moy=np.zeros((nlignes,Tab.shape[1]),dtype=float)
    for j in range(nlignes) :
        moy[j,:]=np.mean(Tab[n*j:(j+1)*n,:],0)
    return moy
    
def moymob(Tab,n) :
    """Fonction faisant une moyenne mobile sur les colonnes d'un tableau
    B=array([[1, 1],
             [2, 2],
             [3, 3],
             [4, 1],
             [4, 1],
             [6, 1]])
    moymob(B,4)=
    array([[ 2.5 ,  1.75],
           [ 3.25,  1.75],
           [ 4.25,  1.5 ]])
    """
    nmoy=np.size(Tab,0)-n+1
    moy=np.zeros((np.size(Tab,1),nmoy),dtype=float)
    for j in range(np.size(Tab,1)) :
        Col=Tab[:,j]
        for i in range(nmoy) :
            moy[j,i]=np.mean(Col[i:i+n])
    return np.transpose(moy)

def rolling_window(a, window):
    """
    Make an ndarray with a rolling window of the last dimension

    Parameters
    ----------
    a : array_like
        Array to add rolling window to
    window : int
        Size of rolling window

    Returns
    -------
    Array that is a view of the original array with a added dimension
    of size w.

    Examples
    --------
    >>> x=np.arange(10).reshape((2,5))
    >>> rolling_window(x, 3)
    array([[[0, 1, 2], [1, 2, 3], [2, 3, 4]],
           [[5, 6, 7], [6, 7, 8], [7, 8, 9]]])

    Calculate rolling mean of last dimension:
    >>> np.mean(rolling_window(x, 3), -1)
    array([[ 1.,  2.,  3.],
           [ 6.,  7.,  8.]])

    """
    if window < 1:
        raise ValueError("`window` must be at least 1.")
    if window > a.shape[-1]:
        raise ValueError("`window` is too long.")
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def moymob2(Tab,n) :
    """Fonction faisant une moyenne mobile sur les colonnes d'un tableau
    B=array([[1, 1],
             [2, 2],
             [3, 3],
             [4, 1],
             [4, 1],
             [6, 1]])
    moymob2(B,4)=
    array([[ 2.5 ,  1.75],
           [ 3.25,  1.75],
           [ 4.25,  1.5 ]])
    """
    return np.transpose(np.mean(rolling_window(Tab.swapaxes(0,1),n),-1))
    
def moymobfull(Tab,n) :
    if n%2==0 :
        p=n/2
    else :
        p=(n-1)/2
    return np.vstack([Tab[:p,],moymob2(Tab,n),Tab[-p:,]])

def movavg1d(Tab,n) :
    if (n-1)%2 != 0 :
        print("nombre invalide pas de la forme 2n+1")
        n=n+1
        print("on prend n=",n)
    p=(n-1)/2
    #~ extended_data = np.hstack([[Tab[0]] * (n- 1), Tab])
    weightings = np.empty((n,), dtype=np.float_)
    weightings[:] = 1.0/n
    #~ np.repeat(1.0, n) / n
    Tabmoy=np.convolve(Tab, weightings, mode='valid')
    Tabstart=np.cumsum(Tab[:p])/np.arange(1,p+1,dtype='float')
    Tabend=np.cumsum(Tab[-p:][::-1])[::-1]/np.arange(p,0,-1,dtype='float')
    return np.r_[Tabstart,Tabmoy,Tabend]
    
def movavg2d(Tab,n) :
    out=np.zeros_like(Tab.T)
    for index,el in enumerate(Tab.T) :
        out[index,:]=movavg1d(el,n)
    
    return out.T

def Pairup(*args) :
    """Fonction permettant de caculer et de stocker l'ensemble des cas
    possibles. Le nombre est laissé libre mais il doit s'agir de 
    tableaux unidimensionnels. Exemple :
    V1=array([1,2])
    V2=array([3,4,5])
    Pairup(V1,V2)=
     array([[1, 3],
            [1, 4],
            [1, 5],
            [2, 3],
            [2, 4],
            [2, 5]])
    """
    M=args[0]
    for j in range(1,len(args)) :
        nM=np.size(M,0)
        nargs=np.size(args[j],0)
        M=np.c_[M.repeat(nargs,axis=0),np.tile(args[j],nM)]
    return M
