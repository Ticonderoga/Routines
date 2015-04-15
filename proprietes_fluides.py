#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math,divers
import numpy as np
try: import psyco; psyco.full() 
except: pass
# Sous programmes de déterminations de propriétés physiques de fluides.
# On se limite ici aux corps simples c.a.d :
# - Air sec
# - Vapeur

#    _____________
#___/ Conversions \____________________________________________________________

def C2K(T):
    """ Conversion d'une température en °C en une température en K."""
    return float(T+273.15)

def K2C(T):
    """ Conversion d'une température en K en une température en °C. """
    if (T<=0) :
        print "K2C (Avertissement) : La température envoyée est négative,\n\
        elle ne peut pas être exprimée en K."
        return float(T)
    else:
        return float(T-273.15)


#    _______________
#___/ hydrogène sec \________________________________________________________________

# Capacité thermique.
def Cph2(T):
    """ Capacité thermique massique (J/g*K) à pression constante de l'h2 sec considéré
    comme un gaz parfait => Cp ne dépend que de la températures "T" exprimée
    en K et comprise entre -50°C et 126°C. """
    # Test de compatibilité de la valeur de la température
    if ((T<223.15) or (T>400)):
        # Message d'erreur et arrêt.
        print " Cph2 (Erreur) : La température envoyée doit être comprise entre\n\
        223.15 K et 400 K."
        return
    else:
        Cpb=[[223.15,224.03,224.91,225.79,226.67,227.55,228.43,229.31,
        230.19,231.07,231.95,232.83,233.71,234.59,235.47,236.35,237.23,
        238.11,238.99,239.87,240.75,241.63,242.51,243.39,244.27,245.15,
        246.03,246.91,247.79,248.67,249.55,250.43,251.31,252.19,253.07,
        253.95,254.83,255.71,256.59,257.47,258.35,259.23,260.11,260.99,
        261.87,262.75,263.63,264.51,265.39,266.27,267.15,268.03,268.91,
        269.79,270.67,271.55,272.43,273.31,274.19,275.07,275.95,276.83,
        277.71,278.59,279.47,280.35,281.23,282.11,282.99,283.87,284.75,
        285.63,286.51,287.39,288.27,289.15,290.03,290.91,291.79,292.67,
        293.55,294.43,295.31,296.19,297.07,297.95,298.83,299.71,300.59,
        301.47,302.35,303.23,304.11,304.99,305.87,306.75,307.63,308.51,
        309.39,310.27,311.15,312.03,312.91,313.79,314.67,315.55,316.43,
        317.31,318.19,319.07,319.95,320.83,321.71,322.59,323.47,324.35,
        325.23,326.11,326.99,327.87,328.75,329.63,330.51,331.39,332.27,
        333.15,334.03,334.91,335.79,336.67,337.55,338.43,339.31,340.19,
        341.07,341.95,342.83,343.71,344.59,345.47,346.35,347.23,348.11,
        348.99,349.87,350.75,351.63,352.51,353.39,354.27,355.15,356.03,
        356.91,357.79,358.67,359.55,360.43,361.31,362.19,363.07,363.95,
        364.83,365.71,366.59,367.47,368.35,369.23,370.11,370.99,371.87,
        372.75,373.63,374.51,375.39,376.27,377.15,378.03,378.91,379.79,
        380.67,381.55,382.43,383.31,384.19,385.07,385.95,386.83,387.71,
        388.59,389.47,390.35,391.23,392.11,392.99,393.87,394.75,395.63,
        396.51,397.39,398.27,399.15
        ],[13.815,13.824,13.833,13.842,13.851,13.86,13.868,13.877,13.885,
        13.893,13.901,13.909,13.917,13.925,13.933,13.941,13.948,13.956,
        13.963,13.971,13.978,13.985,13.992,13.999,14.006,14.013,14.02,
        14.027,14.033,14.04,14.046,14.053,14.059,14.065,14.072,14.078,
        14.084,14.09,14.096,14.101,14.107,14.113,14.118,14.124,14.13,
        14.135,14.14,14.146,14.151,14.156,14.161,14.166,14.171,14.176,
        14.181,14.186,14.19,14.195,14.2,14.204,14.209,14.213,14.218,
        14.222,14.226,14.23,14.235,14.239,14.243,14.247,14.251,14.255,
        14.259,14.262,14.266,14.27,14.274,14.277,14.281,14.284,14.288,
        14.291,14.295,14.298,14.301,14.304,14.308,14.311,14.314,14.317,
        14.32,14.323,14.326,14.329,14.332,14.335,14.337,14.34,14.343,
        14.345,14.348,14.351,14.353,14.356,14.358,14.361,14.363,14.365,
        14.368,14.37,14.372,14.374,14.377,14.379,14.381,14.383,14.385,
        14.387,14.389,14.391,14.393,14.395,14.397,14.399,14.4,14.402,
        14.404,14.406,14.407,14.409,14.411,14.412,14.414,14.415,14.417,
        14.418,14.42,14.421,14.423,14.424,14.426,14.427,14.428,14.429,
        14.431,14.432,14.433,14.434,14.436,14.437,14.438,14.439,14.44,
        14.441,14.442,14.443,14.444,14.445,14.446,14.447,14.448,14.449,
        14.45,14.451,14.452,14.453,14.454,14.454,14.455,14.456,14.457,
        14.458,14.458,14.459,14.46,14.461,14.461,14.462,14.463,14.463,
        14.464,14.464,14.465,14.466,14.466,14.467,14.467,14.468,14.469,
        14.469,14.47,14.47,14.471,14.471,14.472,14.472,14.473,14.473,
        14.473,14.474,14.474]]    
        # 
        i=int(math.floor((T-Cpb[0][0])/(Cpb[0][1]-Cpb[0][0])))
        # On crée les listes de valeurs de températures et de Cp entourant la
        # valeur recherchée.
        LT=[float(0)]*4
        LCp=[float(0)]*4
        for j in range(4):
            LT[j]=Cpb[0][i+j-1]
            LCp[j]=Cpb[1][i+j-1]
        # Puis on appelle la fonction d'interpolation sur ces 4 points
        return divers.interp_polyn_Lagrange(LT,LCp,T)*1000

# Viscosité dynamique.
def muh2(T):
    """ Viscosité dynamique (Pa*s) de l'h2 sec pour une température "T" expimée en K
    et comprise entre -50°C et 126°C."""
    # Test de compatibilité de la valeur de la température
    if ((T<223.15) or (T>400)):
        # Message d'erreur et arrêt.
        print " muas (Erreur) : La température envoyée doit être comprise entre\n\
        223.15 K et 400 K."
        return
    else:
        mub=[[223.15,224.03,224.91,225.79,226.67,227.55,228.43,229.31,
        230.19,231.07,231.95,232.83,233.71,234.59,235.47,236.35,237.23,
        238.11,238.99,239.87,240.75,241.63,242.51,243.39,244.27,245.15,
        246.03,246.91,247.79,248.67,249.55,250.43,251.31,252.19,253.07,
        253.95,254.83,255.71,256.59,257.47,258.35,259.23,260.11,260.99,
        261.87,262.75,263.63,264.51,265.39,266.27,267.15,268.03,268.91,
        269.79,270.67,271.55,272.43,273.31,274.19,275.07,275.95,276.83,
        277.71,278.59,279.47,280.35,281.23,282.11,282.99,283.87,284.75,
        285.63,286.51,287.39,288.27,289.15,290.03,290.91,291.79,292.67,
        293.55,294.43,295.31,296.19,297.07,297.95,298.83,299.71,300.59,
        301.47,302.35,303.23,304.11,304.99,305.87,306.75,307.63,308.51,
        309.39,310.27,311.15,312.03,312.91,313.79,314.67,315.55,316.43,
        317.31,318.19,319.07,319.95,320.83,321.71,322.59,323.47,324.35,
        325.23,326.11,326.99,327.87,328.75,329.63,330.51,331.39,332.27,
        333.15,334.03,334.91,335.79,336.67,337.55,338.43,339.31,340.19,
        341.07,341.95,342.83,343.71,344.59,345.47,346.35,347.23,348.11,
        348.99,349.87,350.75,351.63,352.51,353.39,354.27,355.15,356.03,
        356.91,357.79,358.67,359.55,360.43,361.31,362.19,363.07,363.95,
        364.83,365.71,366.59,367.47,368.35,369.23,370.11,370.99,371.87,
        372.75,373.63,374.51,375.39,376.27,377.15,378.03,378.91,379.79,
        380.67,381.55,382.43,383.31,384.19,385.07,385.95,386.83,387.71,
        388.59,389.47,390.35,391.23,392.11,392.99,393.87,394.75,395.63,
        396.51,397.39,398.27,399.15],
        [7.3098E-06,7.3296E-06,7.3493E-06,7.3691E-06,7.3888E-06,
        7.4085E-06,7.4281E-06,7.4478E-06,7.4674E-06,7.487E-06,7.5065E-06,
        7.5261E-06,7.5456E-06,7.5651E-06,7.5845E-06,7.604E-06,7.6234E-06,
        7.6428E-06,7.6621E-06,7.6815E-06,7.7008E-06,7.7201E-06,7.7394E-06,
        7.7586E-06,7.7779E-06,7.7971E-06,7.8163E-06,7.8354E-06,7.8546E-06,
        7.8737E-06,7.8928E-06,7.9119E-06,7.9309E-06,7.9499E-06,7.969E-06,
        7.9879E-06,8.0069E-06,8.0258E-06,8.0448E-06,8.0637E-06,8.0825E-06,
        8.1014E-06,8.1202E-06,8.1391E-06,8.1578E-06,8.1766E-06,8.1954E-06,
        8.2141E-06,8.2328E-06,8.2515E-06,8.2702E-06,8.2888E-06,8.3074E-06,
        8.3261E-06,8.3446E-06,8.3632E-06,8.3818E-06,8.4003E-06,8.4188E-06,
        8.4373E-06,8.4558E-06,8.4742E-06,8.4926E-06,8.511E-06,8.5294E-06,
        8.5478E-06,8.5662E-06,8.5845E-06,8.6028E-06,8.6211E-06,8.6394E-06,
        8.6576E-06,8.6759E-06,8.6941E-06,8.7123E-06,8.7305E-06,8.7486E-06,
        8.7668E-06,8.7849E-06,8.803E-06,8.8211E-06,8.8392E-06,8.8572E-06,
        8.8753E-06,8.8933E-06,8.9113E-06,8.9292E-06,8.9472E-06,8.9652E-06,
        8.9831E-06,9.001E-06,9.0189E-06,9.0368E-06,9.0546E-06,9.0725E-06,
        9.0903E-06,9.1081E-06,9.1259E-06,9.1436E-06,9.1614E-06,9.1791E-06,
        9.1969E-06,9.2146E-06,9.2322E-06,9.2499E-06,9.2676E-06,9.2852E-06,
        9.3028E-06,9.3204E-06,9.338E-06,9.3556E-06,9.3731E-06,9.3907E-06,
        9.4082E-06,9.4257E-06,9.4432E-06,9.4607E-06,9.4781E-06,9.4956E-06,
        9.513E-06,9.5304E-06,9.5478E-06,9.5652E-06,9.5825E-06,9.5999E-06,
        9.6172E-06,9.6345E-06,9.6518E-06,9.6691E-06,9.6864E-06,9.7036E-06,
        9.7209E-06,9.7381E-06,9.7553E-06,9.7725E-06,9.7897E-06,9.8068E-06,
        9.824E-06,9.8411E-06,9.8582E-06,9.8753E-06,9.8924E-06,9.9095E-06,
        9.9265E-06,9.9436E-06,9.9606E-06,9.9776E-06,9.9946E-06,1.0012E-05,
        1.0029E-05,1.0046E-05,1.0062E-05,1.0079E-05,1.0096E-05,1.0113E-05,
        1.013E-05,1.0147E-05,1.0164E-05,1.0181E-05,1.0198E-05,1.0214E-05,
        1.0231E-05,1.0248E-05,1.0265E-05,1.0281E-05,1.0298E-05,1.0315E-05,
        1.0332E-05,1.0348E-05,1.0365E-05,1.0382E-05,1.0398E-05,1.0415E-05,
        1.0432E-05,1.0448E-05,1.0465E-05,1.0482E-05,1.0498E-05,1.0515E-05,
        1.0531E-05,1.0548E-05,1.0564E-05,1.0581E-05,1.0597E-05,1.0614E-05,
        1.063E-05,1.0647E-05,1.0663E-05,1.068E-05,1.0696E-05,1.0713E-05,
        1.0729E-05,1.0745E-05,1.0762E-05,1.0778E-05,1.0794E-05,1.0811E-05,
        1.0827E-05,1.0843E-05,1.086E-05,1.0876E-05]]
        #
        i=int(math.floor((T-mub[0][0])/(mub[0][1]-mub[0][0])))
        # On crée les listes de valeurs de températures et de Cp entourant la
        # valeur recherchée.
        LT=[float(0)]*4
        Lmu=[float(0)]*4
        for j in range(4):
            LT[j]=mub[0][i+j-1]
            Lmu[j]=mub[1][i+j-1]
        # Puis on appelle la fonction d'interpolation sur ces 4 points
        return divers.interp_polyn_Lagrange(LT,Lmu,T)

# Viscosité cinématique.
def nuh2(p,T):
    """ Viscosité cinématique de l'h2 sec pour une pression "p" donnée et une
    température "T" exprimée en K et comprise entre -50°C et 126°C. """
    # Test de compatibilité de la valeur de la température
    if ((T<223.15) or (T>400)):
        # Message d'erreur et arrêt.
        print " nuas (Erreur) : La température envoyée doit être comprise entre\n\
        223.15 K et 400 K."
        return
    else:
        # On commence par calculer la valeur de la viscosité dynamique
        mu=muh2(T)
        # Puis on calcule la masse volumique par la loi d'état des gaz parfaits, ce
        # qui nous donne la valeur de la viscosité cinématique.
        return mu/(p/(4124.*T))

# Conductivité thermique.
def lambdah2(T):
    """ Conductivité thermique de l'h2 sec pour une température "T" exprimée
    en K et comprise entre -50°C et 126°C."""
    # Test de compatibilité de la valeur de la température
    if ((T<223.15) or (T>400)):
        # Message d'erreur et arrêt.
        print " lambdaas (Erreur) : La température envoyée doit être comprise entre\n\
        223.15 K et 573.15 K."
        return
    else:
        lambdab=[[223.15,224.03,224.91,225.79,226.67,227.55,228.43,229.31,
        230.19,231.07,231.95,232.83,233.71,234.59,235.47,236.35,237.23,
        238.11,238.99,239.87,240.75,241.63,242.51,243.39,244.27,245.15,
        246.03,246.91,247.79,248.67,249.55,250.43,251.31,252.19,253.07,
        253.95,254.83,255.71,256.59,257.47,258.35,259.23,260.11,260.99,
        261.87,262.75,263.63,264.51,265.39,266.27,267.15,268.03,268.91,
        269.79,270.67,271.55,272.43,273.31,274.19,275.07,275.95,276.83,
        277.71,278.59,279.47,280.35,281.23,282.11,282.99,283.87,284.75,
        285.63,286.51,287.39,288.27,289.15,290.03,290.91,291.79,292.67,
        293.55,294.43,295.31,296.19,297.07,297.95,298.83,299.71,300.59,
        301.47,302.35,303.23,304.11,304.99,305.87,306.75,307.63,308.51,
        309.39,310.27,311.15,312.03,312.91,313.79,314.67,315.55,316.43,
        317.31,318.19,319.07,319.95,320.83,321.71,322.59,323.47,324.35,
        325.23,326.11,326.99,327.87,328.75,329.63,330.51,331.39,332.27,
        333.15,334.03,334.91,335.79,336.67,337.55,338.43,339.31,340.19,
        341.07,341.95,342.83,343.71,344.59,345.47,346.35,347.23,348.11,
        348.99,349.87,350.75,351.63,352.51,353.39,354.27,355.15,356.03,
        356.91,357.79,358.67,359.55,360.43,361.31,362.19,363.07,363.95,
        364.83,365.71,366.59,367.47,368.35,369.23,370.11,370.99,371.87,
        372.75,373.63,374.51,375.39,376.27,377.15,378.03,378.91,379.79,
        380.67,381.55,382.43,383.31,384.19,385.07,385.95,386.83,387.71,
        388.59,389.47,390.35,391.23,392.11,392.99,393.87,394.75,395.63,
        396.51,397.39,398.27,399.15
        ],[0.14584,0.14633,0.14683,0.14733,0.14782,0.14832,0.14881,
        0.14931,0.1498,0.15029,0.15078,0.15126,0.15175,0.15224,0.15273,
        0.15321,0.1537,0.15419,0.15468,0.15516,0.15564,0.15611,0.15658,
        0.15705,0.15752,0.15799,0.15846,0.15893,0.1594,0.15987,0.16035,
        0.16081,0.16128,0.16174,0.1622,0.16266,0.16313,0.16359,0.16405,
        0.16451,0.16498,0.16544,0.1659,0.16636,0.16681,0.16727,0.16772,
        0.16817,0.16863,0.16908,0.16954,0.16999,0.17045,0.1709,0.17135,
        0.17178,0.17222,0.17266,0.1731,0.17353,0.17397,0.17441,0.17485,
        0.17529,0.17572,0.17616,0.1766,0.17704,0.17748,0.17792,0.17835,
        0.17879,0.17923,0.17967,0.18011,0.18055,0.18099,0.18141,0.18183,
        0.18225,0.18267,0.18309,0.18352,0.18394,0.18436,0.18478,0.1852,
        0.18563,0.18605,0.18647,0.18689,0.18731,0.18774,0.18816,0.18858,
        0.189,0.18943,0.18985,0.19027,0.19069,0.19109,0.1915,0.1919,
        0.19231,0.19271,0.19312,0.19352,0.19393,0.19434,0.19474,0.19515,
        0.19559,0.19604,0.19649,0.19694,0.19739,0.19784,0.19829,0.19874,
        0.19919,0.19964,0.20009,0.20054,0.20098,0.20142,0.20186,0.2023,
        0.20274,0.20319,0.20363,0.20407,0.20451,0.20495,0.20539,0.20584,
        0.20628,0.20672,0.20716,0.2076,0.20805,0.20849,0.20893,0.20937,
        0.20982,0.21026,0.21069,0.21113,0.21156,0.212,0.21243,0.21286,
        0.2133,0.21373,0.21417,0.2146,0.21503,0.21546,0.21586,0.21627,
        0.21668,0.21709,0.2175,0.21791,0.21831,0.21872,0.21913,0.21954,
        0.21995,0.22038,0.2208,0.22123,0.22166,0.22208,0.22251,0.22293,
        0.22336,0.22379,0.22421,0.22464,0.22506,0.22546,0.22587,0.22628,
        0.22669,0.2271,0.22751,0.22792,0.22833,0.22874,0.22915,0.22956,
        0.22997,0.23038,0.23079,0.2312,0.23161,0.23202,0.23243,0.23284,
        0.23325,0.23366,]]
        # 
        i=int(math.floor((T-lambdab[0][0])/(lambdab[0][1]-lambdab[0][0])))
        # On crée les listes de valeurs de températures et de lambda entourant la
        # valeur recherchée.
        LT=[float(0)]*4
        Ll=[float(0)]*4
        for j in range(4):
            LT[j]=lambdab[0][i+j-1]
            Ll[j]=lambdab[1][i+j-1]
        # Puis on appelle la fonction d'interpolation sur ces 4 points
        
        return divers.interp_polyn_Lagrange(LT,Ll,T)

# masse volumique.
def rhoh2(T):
    """ masse volumique de l'h2 sec pour une température "T" exprimée
    en K et comprise entre -50°C et 126°C."""
    # Test de compatibilité de la valeur de la température
    if ((T<223.15) or (T>399.15)):
        # Message d'erreur et arrêt.
        print " rhoh2 (Erreur) : La température envoyée doit être comprise entre\n\
        223.15 K et 400 K."
        return
    else:
        rhob=[[223.15,224.03,224.91,225.79,226.67,227.55,228.43,229.31,
        230.19,231.07,231.95,232.83,233.71,234.59,235.47,236.35,237.23,
        238.11,238.99,239.87,240.75,241.63,242.51,243.39,244.27,245.15,
        246.03,246.91,247.79,248.67,249.55,250.43,251.31,252.19,253.07,
        253.95,254.83,255.71,256.59,257.47,258.35,259.23,260.11,260.99,
        261.87,262.75,263.63,264.51,265.39,266.27,267.15,268.03,268.91,
        269.79,270.67,271.55,272.43,273.31,274.19,275.07,275.95,276.83,
        277.71,278.59,279.47,280.35,281.23,282.11,282.99,283.87,284.75,
        285.63,286.51,287.39,288.27,289.15,290.03,290.91,291.79,292.67,
        293.55,294.43,295.31,296.19,297.07,297.95,298.83,299.71,300.59,
        301.47,302.35,303.23,304.11,304.99,305.87,306.75,307.63,308.51,
        309.39,310.27,311.15,312.03,312.91,313.79,314.67,315.55,316.43,
        317.31,318.19,319.07,319.95,320.83,321.71,322.59,323.47,324.35,
        325.23,326.11,326.99,327.87,328.75,329.63,330.51,331.39,332.27,
        333.15,334.03,334.91,335.79,336.67,337.55,338.43,339.31,340.19,
        341.07,341.95,342.83,343.71,344.59,345.47,346.35,347.23,348.11,
        348.99,349.87,350.75,351.63,352.51,353.39,354.27,355.15,356.03,
        356.91,357.79,358.67,359.55,360.43,361.31,362.19,363.07,363.95,
        364.83,365.71,366.59,367.47,368.35,369.23,370.11,370.99,371.87,
        372.75,373.63,374.51,375.39,376.27,377.15,378.03,378.91,379.79,
        380.67,381.55,382.43,383.31,384.19,385.07,385.95,386.83,387.71,
        388.59,389.47,390.35,391.23,392.11,392.99,393.87,394.75,395.63,
        396.51,397.39,398.27,399.15
        ],[0.11002,0.10959,0.10916,0.10873,0.10831,0.10789,0.10747,
        0.10706,0.10665,0.10625,0.10584,0.10544,0.10505,0.10465,0.10426,
        0.10387,0.10349,0.10311,0.10273,0.10235,0.10198,0.10161,0.10124,
        0.10087,0.10051,0.10015,0.099789,0.099433,0.09908,0.09873,
        0.098381,0.098036,0.097693,0.097352,0.097013,0.096677,0.096344,
        0.096012,0.095683,0.095356,0.095031,0.094709,0.094389,0.09407,
        0.093754,0.093441,0.093129,0.092819,0.092511,0.092206,0.091902,
        0.091601,0.091301,0.091003,0.090707,0.090414,0.090122,0.089832,
        0.089543,0.089257,0.088972,0.08869,0.088409,0.08813,0.087852,
        0.087577,0.087303,0.08703,0.08676,0.086491,0.086224,0.085958,
        0.085694,0.085432,0.085171,0.084912,0.084655,0.084399,0.084144,
        0.083891,0.08364,0.08339,0.083142,0.082895,0.082649,0.082405,
        0.082163,0.081922,0.081682,0.081444,0.081207,0.080971,0.080737,
        0.080504,0.080273,0.080042,0.079813,0.079586,0.07936,0.079135,
        0.078911,0.078688,0.078467,0.078247,0.078029,0.077811,0.077595,
        0.07738,0.077166,0.076953,0.076741,0.076531,0.076322,0.076114,
        0.075907,0.075701,0.075496,0.075292,0.07509,0.074888,0.074688,
        0.074489,0.074291,0.074093,0.073897,0.073702,0.073508,0.073315,
        0.073123,0.072932,0.072742,0.072553,0.072365,0.072177,0.071991,
        0.071806,0.071622,0.071439,0.071256,0.071075,0.070894,0.070715,
        0.070536,0.070358,0.070181,0.070005,0.06983,0.069656,0.069483,
        0.06931,0.069138,0.068968,0.068798,0.068628,0.06846,0.068293,
        0.068126,0.06796,0.067795,0.067631,0.067467,0.067305,0.067143,
        0.066982,0.066821,0.066662,0.066503,0.066345,0.066188,0.066031,
        0.065875,0.06572,0.065566,0.065412,0.065259,0.065107,0.064956,
        0.064805,0.064655,0.064505,0.064357,0.064209,0.064061,0.063915,
        0.063769,0.063623,0.063479,0.063335,0.063191,0.063048,0.062906,
        0.062765,0.062624,0.062484,0.062344,0.062206,0.062067,0.06193,
        0.061792,0.061656,0.06152
        ]]
        # 
        i=int(math.floor((T-rhob[0][0])/(rhob[0][0]-rhob[0][1])))
        # On crée les listes de valeurs de températures et de lambda entourant la
        # valeur recherchée.
        LT=[float(0)]*4
        Ll=[float(0)]*4
        for j in range(4):
            LT[j]=rhob[0][i+j-1]
            Ll[j]=rhob[1][i+j-1]
        # Puis on appelle la fonction d'interpolation sur ces 4 points
        return divers.interp_polyn_Lagrange(LT,Ll,T)


#    _________
#___/ Air sec \________________________________________________________________

# Capacité thermique.
def Cpas(T):
    """ Capacité thermique massique à pression constante de l'air sec considéré
    comme un gaz parfait => Cp ne dépend que de la températures "T" exprimée
    en K et comprise entre -50°C et 300°C. """
    # Test de compatibilité de la valeur de la température
    if ((T<223.15) or (T>573.15)):
        # Message d'erreur et arrêt.
        print " Cpas (Erreur) : La température envoyée doit être comprise entre\n\
        223.15 K et 573.15 K."
        return
    else:
        Cpb=[[
        223.15,233.15,243.15,253.15,263.15,273.15,283.15,293.15,303.15,313.15,
        323.15,333.15,343.15,353.15,363.15,373.15,383.15,393.15,403.15,413.15,
        423.15,433.15,443.15,453.15,463.15,473.15,483.15,493.15,503.15,513.15,
        523.15,533.15,543.15,553.15,563.15,573.15
        ],[
        1006.4,1006,1005.8,1005.7,1005.6,1005.7,1005.8,1006.1,1006.4,1006.8,
        1007.4,1008,1008.7,1009.5,1010.3,1011.3,1012.3,1013.4,1014.6,1015.9,
        1017.2,1018.6,1020.1,1021.7,1023.3,1025,1026.8,1028.6,1030.5,1032.4,
        1034.4,1036.5,1038.6,1040.7,1042.9,1045.2]]
        # 
        i=int(math.floor(float(T-Cpb[0][0])/10))
        # On crée les listes de valeurs de températures et de Cp entourant la
        # valeur recherchée.
        LT=[float(0)]*4
        LCp=[float(0)]*4
        for j in range(4):
            LT[j]=Cpb[0][i+j-1]
            LCp[j]=Cpb[1][i+j-1]
        # Puis on appelle la fonction d'interpolation sur ces 4 points
        return divers.interp_polyn_Lagrange(LT,LCp,T)

# Viscosité dynamique.
def muas(T):
    """ Viscosité dynamique de l'air sec pour une température "T" expimée en K
    et comprise entre -50°C et 300°C."""
    # Test de compatibilité de la valeur de la température
    if ((T<223.15) or (T>573.15)):
        # Message d'erreur et arrêt.
        print " muas (Erreur) : La température envoyée doit être comprise entre\n\
        223.15 K et 573.15 K."
        return
    else:
        mub=[[223.15,233.15,243.15,253.15,263.15,273.15,283.15,293.15,303.15,
        313.15,323.15,333.15,343.15,353.15,363.15,373.15,383.15,393.15,403.15,
        413.15,423.15,433.15,443.15,453.15,463.15,473.15,483.15,493.15,503.15,
        513.15,523.15,533.15,543.15,553.15,563.15,573.15],[14.63e-6,15.17e-6,
        15.69e-6,16.20e-6,16.71e-6,17.20e-6,17.69e-6,18.17e-6,18.65e-6,
        19.11e-6,19.57e-6,20.03e-6,20.47e-6,20.92e-6,21.35e-6,21.78e-6,
        22.20e-6,22.62e-6,23.03e-6,23.44e-6,23.84e-6,24.24e-6,24.63e-6,
        25.03e-6,25.41e-6,25.79e-6,26.17e-6,26.54e-6,26.91e-6,27.27e-6,
        27.64e-6,27.99e-6,28.35e-6,28.70e-6,29.05e-6,29.39e-6]]
        #
        i=int(math.floor(float(T-mub[0][0])/10))
        # On crée les listes de valeurs de températures et de Cp entourant la
        # valeur recherchée.
        LT=[float(0)]*4
        Lmu=[float(0)]*4
        for j in range(4):
            LT[j]=mub[0][i+j-1]
            Lmu[j]=mub[1][i+j-1]
        # Puis on appelle la fonction d'interpolation sur ces 4 points
        return divers.interp_polyn_Lagrange(LT,Lmu,T)

# Viscosité cinématique.
def nuas(p,T):
    """ Viscosité cinématique de l'air sec pour une pression "p" donnée et une
    température "T" exprimée en K et comprise entre -50°C et 300°C. """
    # Test de compatibilité de la valeur de la température
    if ((T<223.15) or (T>573.15)):
        # Message d'erreur et arrêt.
        print " nuas (Erreur) : La température envoyée doit être comprise entre\n\
        223.15 K et 573.15 K."
        return
    else:
        # On commence par calculer la valeur de la viscosité dynamique
        mu=muas(T)
        # Puis on calcule la masse volumique par la loi d'état des gaz parfaits, ce
        # qui nous donne la valeur de la viscosité cinématique.
        return mu/(p/(287.*T))

# Conductivité thermique.
def lambdaas(T):
    """ Conductivité thermique de l'air sec pour une température "T" exprimée
    en K et comprise entre -50°C et 300°C."""
    # Test de compatibilité de la valeur de la température
    if ((T<223.15) or (T>573.15)):
        # Message d'erreur et arrêt.
        print " lambdaas (Erreur) : La température envoyée doit être comprise entre\n\
        223.15 K et 573.15 K."
        return
    else:
        lambdab=[[
        223.15,233.15,243.15,253.15,263.15,273.15,283.15,293.15,303.15,313.15,
        323.15,333.15,343.15,353.15,363.15,373.15,383.15,393.15,403.15,413.15,
        423.15,433.15,443.15,453.15,463.15,473.15,483.15,493.15,503.15,513.15,
        523.15,533.15,543.15,553.15,563.15,573.15
        ],[
        20.04e-3,20.86e-3,21.68e-3,22.49e-3,23.29e-3,24.08e-3,24.87e-3,
        25.64e-3,26.38e-3,27.10e-3,27.81e-3,28.52e-3,29.22e-3,29.91e-3,
        30.59e-3,31.27e-3,31.94e-3,32.61e-3,33.28e-3,33.94e-3,34.59e-3,
        35.25e-3,35.89e-3,36.54e-3,37.18e-3,37.81e-3,38.45e-3,39.08e-3,
        39.71e-3,40.33e-3,40.95e-3,41.57e-3,42.18e-3,42.79e-3,43.40e-3,
        44.01e-3]]
        # 
        i=int(math.floor(float(T-lambdab[0][0])/10))
        # On crée les listes de valeurs de températures et de lambda entourant la
        # valeur recherchée.
        LT=[float(0)]*4
        Ll=[float(0)]*4
        for j in range(4):
            LT[j]=lambdab[0][i+j-1]
            Ll[j]=lambdab[1][i+j-1]
        # Puis on appelle la fonction d'interpolation sur ces 4 points
        return divers.interp_polyn_Lagrange(LT,Ll,T)

# Nombre de Prandl
def Pras(T):
    """ Nombre de Prandtl de l'air sec pour une température "T" exprimée en K
    et comprise entre -50 et 300 °C """
    # Test de compatibilité de la valeur de la température
    if ((T<223.15) or (T>573.15)):
        # Message d'erreur et arrêt.
        print " Pras (Erreur) : La température envoyée doit être comprise entre\n\
        223.15 K et 573.15 K."
        return
    else:
        Prb=[[
        223.15,233.15,243.15,253.15,263.15,273.15,283.15,293.15,303.15,313.15,
        323.15,333.15,343.15,353.15,363.15,373.15,383.15,393.15,403.15,413.15,
        423.15,433.15,443.15,453.15,463.15,473.15,483.15,493.15,503.15,513.15,
        523.15,533.15,543.15,553.15,563.15,573.15
        ],[
        0.735,0.731,0.728,0.724,0.721,0.718,0.716,0.713,0.712,0.710,0.709,
        0.708,0.707,0.706,0.705,0.704,0.704,0.703,0.702,0.702,0.701,0.701,
        0.700,0.700,0.699,0.699,0.699,0.699,0.698,0.698,0.698,0.698,0.698,
        0.698,0.698,0.698]]
        # 
        i=int(math.floor(float(T-Prb[0][0])/10))
        # On crée les listes de valeurs de températures et de lambda entourant la
        # valeur recherchée.
        LT=[float(0)]*4
        LPr=[float(0)]*4
        for j in range(4):
            LT[j]=Prb[0][i+j-1]
            LPr[j]=Prb[1][i+j-1]
        # Puis on appelle la fonction d'interpolation sur ces 4 points
        return divers.interp_polyn_Lagrange(LT,LPr,T)

#    ______________
#___/ Vapeur d'eau \___________________________________________________________

# Capacité thermique.
def Cpvap(T):
    """ Capacité thermique massique à pression constante de la vapeur d'eau 
    considéréé comme un gaz parfait => Cp ne dépend que de la températures "T"
    exprimée en K et comprise entre 0°C et 300°C. """
    # Test de compatibilité de la valeur de la température
    if ((T<273.15) or (T>573.15)):
        # Message d'erreur et arrêt.
        print " Cpvap (Erreur) : La température envoyée doit être comprise entre\n\
        273.15 K et 573.15 K."
        return
    else:
        Cpb=[[
        273.15,283.15,293.15,303.15,313.15,323.15,333.15,343.15,353.15,363.15,
        373.15,383.15,393.15,403.15,413.15,423.15,433.15,443.15,453.15,463.15,
        473.15,483.15,493.15,503.15,513.15,523.15,533.15,543.15,553.15,563.15,
        573.15],[
        1.8844,1.8947,1.9059,1.9180,1.9314,1.9468,1.9648,1.9862,2.0120,2.0430,
        2.0801,2.1244,2.1770,2.2389,2.3110,2.3940,2.4884,2.5946,2.7130,2.8444,
        2.9897,3.1505,3.3291,3.5287,3.7539,4.0108,4.3079,4.6566,5.0736,5.5826,
        6.2204]]
        # 
        i=int(math.floor(float(T-Cpb[0][0])/10))
        # On crée les listes de valeurs de températures et de Cp entourant la
        # valeur recherchée.
        LT=[float(0)]*4
        LCp=[float(0)]*4
        for j in range(4):
            LT[j]=Cpb[0][i+j-1]
            LCp[j]=Cpb[1][i+j-1]
        # Puis on appelle la fonction d'interpolation sur ces 4 points
        return divers.interp_polyn_Lagrange(LT,LCp,T)*1e3

# Viscosité dynamique.
def muvap(T):
    """ Viscosité dynamique de la vapeur d'eau pour une température "T" exprimée
    en K et comprise entre 0°C et 300°C."""
    # Test de compatibilité de la valeur de la température
    if ((T<273.15) or (T>573.15)):
        # Message d'erreur et arrêt.
        print " muvap (Erreur) : La température envoyée doit être comprise entre\n\
        273.15 K et 573.15 K."
        return
    else:
        mub=[[
        273.15,283.15,293.15,303.15,313.15,323.15,333.15,343.15,353.15,363.15,
        373.15,383.15,393.15,403.15,413.15,423.15,433.15,443.15,453.15,463.15,
        473.15,483.15,493.15,503.15,513.15,523.15,533.15,543.15,553.15,563.15,
        573.15],[
        9.2163e-06,9.4614e-06,9.7275e-06,1.0011e-05,1.0308e-05,1.0617e-05,
        1.0935e-05,1.1261e-05,1.1593e-05,1.1929e-05,1.2270e-05,1.2612e-05,
        1.2956e-05,1.3302e-05,1.3647e-05,1.3992e-05,1.4337e-05,1.4682e-05,
        1.5026e-05,1.5370e-05,1.5715e-05,1.6062e-05,1.6411e-05,1.6765e-05,
        1.7125e-05,1.7495e-05,1.7877e-05,1.8277e-05,1.8700e-05,1.9155e-05,
        1.9652e-05]]
        # 
        i=int(math.floor(float(T-mub[0][0])/10))
        # On crée les listes de valeurs de températures et de mu entourant la
        # valeur recherchée.
        LT=[float(0)]*4
        Lmu=[float(0)]*4
        for j in range(4):
            LT[j]=mub[0][i+j-1]
            Lmu[j]=mub[1][i+j-1]
        # Puis on appelle la fonction d'interpolation sur ces 4 points
        return divers.interp_polyn_Lagrange(LT,Lmu,T)

# Viscosité cinématique.
def nuvap(p,T):
    """ Viscosité cinématique de la vapeur d'eau pour une pression "p" donnée
    et une température "T" exprimée en K et comprise entre 0°C et 300°C. """
    # Test de compatibilité de la valeur de la température
    if ((T<273.15) or (T>573.15)):
        # Message d'erreur et arrêt.
        print " nuvap (Erreur) : La température envoyée doit être comprise entre\n\
        273.15 K et 573.15 K."
        return
    else:
        # On commence par calculer la valeur de la viscosité dynamique
        mu=muvap(T)
        # Puis on calcule la masse volumique par la loi d'état des gaz parfaits, ce
        # qui nous donne la valeur de la viscosité cinématique.
        return mu/(p/(462.*T))

# Conductivité thermique.
def lambdavap(T):
    """ Conductivité thermique de la vapeur d'eau pour une température "T"
    exprimée en K et comprise entre 0°C et 300°C."""
    # Test de compatibilité de la valeur de la température
    if ((T<223.15) or (T>573.15)):
        # Message d'erreur et arrêt.
        print " lambdavap (Erreur) : La température envoyée doit être comprise entre\n\
        223.15 K et 573.15 K."
        return
    else:
        lambdab=[[
        223.15,233.15,243.15,253.15,263.15,273.15,283.15,293.15,303.15,313.15,
        323.15,333.15,343.15,353.15,363.15,373.15,383.15,393.15,403.15,413.15,
        423.15,433.15,443.15,453.15,463.15,473.15,483.15,493.15,503.15,513.15,
        523.15,533.15,543.15,553.15,563.15,573.15
        ],[
        0.017071,0.017622,0.018228,0.018887,0.019599,0.020365,0.021188,0.022069,
        0.023012,0.024020,0.025097,0.026246,0.027468,0.028766,0.030142,0.031597,
        0.033132,0.034750,0.036450,0.038238,0.040115,0.042090,0.044172,0.046378,
        0.048730,0.051266,0.054035,0.057115,0.060618,0.064713,0.069655]]
        # 
        i=int(math.floor(float(T-lambdab[0][0])/10))
        # On crée les listes de valeurs de températures et de lambda entourant la
        # valeur recherchée.
        LT=[float(0)]*4
        Ll=[float(0)]*4
        for j in range(4):
            LT[j]=lambdab[0][i+j-1]
            Ll[j]=lambdab[1][i+j-1]
        # Puis on appelle la fonction d'interpolation sur ces 4 points
        return divers.interp_polyn_Lagrange(LT,Ll,T)
        

# Nombre de Prandl
def Prvap(T):
    """ Nombre de Prandtl de la vapeur d'eau pour une température "T" exprimée 
    en K et comprise entre 0 et 300 °C """
    # Test de compatibilité de la valeur de la température
    if ((T<273.15) or (T>573.15)):
        # Message d'erreur et arrêt.
        print " Prvap (Erreur) : La température envoyée doit être comprise entre\n\
        273.15 K et 573.15 K."
        return
    else:
        Prb=[[
        273.15,283.15,293.15,303.15,313.15,323.15,333.15,343.15,353.15,363.15,
        373.15,383.15,393.15,403.15,413.15,423.15,433.15,443.15,453.15,463.15,
        473.15,483.15,493.15,503.15,513.15,523.15,533.15,543.15,553.15,563.15,
        573.15],[
        1.0173,1.0173,1.0171,1.0166,1.0158,1.0149,1.0140,1.0139,1.0136,1.0146,
        1.0170,1.0208,1.0268,1.0353,1.0463,1.0601,1.0768,1.0962,1.1184,1.1433,
        1.1712,1.2023,1.2368,1.2756,1.3192,1.3687,1.4252,1.4901,1.5651,1.6524,
        1.7550]]
        # 
        i=int(math.floor(float(T-Prb[0][0])/10))
        # On crée les listes de valeurs de températures et de lambda entourant la
        # valeur recherchée.
        LT=[float(0)]*4
        LPr=[float(0)]*4
        for j in range(4):
            LT[j]=Prb[0][i+j-1]
            LPr[j]=Prb[1][i+j-1]
        # Puis on appelle la fonction d'interpolation sur ces 4 points
        return divers.interp_polyn_Lagrange(LT,LPr,T)

# Chaleur latente de vaporisation
def ChLv(T):
    """ Chaleur latente de Vaporisation en J/kg avec la température 
    en K et comprise entre 0 et 300 °C """
    # Test de compatibilité de la valeur de la température
    if ((T<273.16) or (T>573.15)):
        # Message d'erreur et arrêt.
        print " Lv (Erreur) : La température envoyée doit être comprise entre\n\
        273.16 K et 573.16 K."
        return
    else:
        Lvb=[[
        273.16,283.16,293.16,303.16,313.16,323.16,333.16,343.16,353.16,
        363.16,373.16,383.16,393.16,403.16,413.16,423.16,433.16,443.16,
        453.16,463.16,473.16,483.16,493.16,503.16,513.16,523.16,533.16,
        543.16,553.16,563.16,573.16],[
        2500.9,2477.14,2453.54,2429.82,2405.93,2381.92,2357.68,2332.99,
        2307.95,2282.42,2256.39,2229.64,2202.05,2173.67,2144.3,2113.68,
        2081.98,2048.77,2014.1,1977.82,1939.68,1899.62,1857.28,1812.66,
        1765.4,1715.1,1661.6,1604.4,1543,1476.6,1404.5]]
        # 
        i=int(math.floor(float(T-Lvb[0][0])/10))
        # On crée les listes de valeurs de températures et de Lv entourant la
        # valeur recherchée.
        LT=[float(0)]*4
        LLv=[float(0)]*4
        for j in range(4):
            LT[j]=Lvb[0][i+j-1]
            LLv[j]=Lvb[1][i+j-1]
        # Puis on appelle la fonction d'interpolation sur ces 4 points
        return divers.interp_polyn_Lagrange(LT,LLv,T)*1e3


# Pression de saturation.
def ps(T):
    """ Calcul de la pression de saturation de vapeur à "T" donnée en K. """
    if (T >= 273.16):
        C=[22105649.25,-27405.526,97.5413,-0.146244,0.12558E-3,-0.48502E-7,4.34903,-0.39381E-2]
        if (T >= 534.):
            print "ps (Avertissement) : La correlation utilisee pour calculer\n\
            la pression de vapeur saturante n'est normalement valable que \n\
            jusqu'a 260 C."
        return float(C[0]*math.exp((C[1]+C[2]*T+C[3]*(T**2)+C[4]*(T**3)+C[5]*(T**4))/(C[6]*T+C[7]*(T**2))))
    else:
        if (T < 243.):
            print "ps (Avertissement) : La correlation utilisee pour calculer\n\
            la pression de vapeur saturante n'est normalement valable qu'a \n\
            partir de -18 C.",T
        T2=K2C(T)
        return float(math.exp((1513.+25.6*T2)/(236.+T2)))
            # Pv=exp((1513+25.6*TC)/(236+TC))
            
def Pvap(TC) :
    """ Calcul de la pression de saturation de vapeur à "T" donnée en C. """
    return ps(C2K(TC))


# Masse volumique
def rhoas(p,T):
    """Calcul de la masse volumique d'un air sec à pression "p" (en Pa) 
    à température "T" (en K) le résultat est kg/m^³"""
    return 1/T*p/287.05

# Coefficient de diffusion de la vapeur d'eau dans l'air.
def Dasvap(p,T):
    """Coefficient de diffusion de masse de la vapeur d'eau dans l'air pour "p"
    exprimée en Pa et "T" en K. -20°C < T < 300°C. """
    if (T < 253.15):
        print "Dasvap (Erreur) : La température doit être comprise entre -20°C\n\
        et 300°C."
        return
    elif (T >= 253.15) and (T < 353.15):
        return 104.91143e-6*(T**(1.774))/float(p)
    elif (T >= 353.15) and (T < 573.15):
        return 805.2375e-6/float(p)*(T**(5./2))/(190+T)
    else:
        print "Dasvap (Erreur) : La température doit être comprise entre -20°C\n\
        et 300°C."
        return 
