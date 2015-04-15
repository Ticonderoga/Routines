# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 15:55:26 2010

@author: -
"""

import DiagHumid as Humid
import matplotlib.pyplot as plt

plt.close('all')
A_bg=Humid.air()
A_bg.definir(["T","HR"],[10.0,0.05])
A_hd=Humid.air()
A_hd.definir(["T","na"],[150.0,0.04])
LH,Lna,LT,LHR,Lair=Humid.Init_diag(A_bg,A_hd,1)
Humid.Trace_diag(LH,Lna,LT,LHR,Lair,1)

#~ AK2SO4=Humid.air()
#~ AK2SO4.definir(["T","Tr"],[19.7,19.0])
#~ Humid.Fill_Air(AK2SO4)
#~ AK2SO4.affiche()

#~ x,y=Humid.Transf_xy(AK2SO4.car[1],AK2SO4.car[0])
#~ plt.plot([x],[y],'ro')
#~ CpK2S04=Humid.Cpah(AK2SO4)
#~ PrK2S04=Humid.Prah(AK2SO4)
