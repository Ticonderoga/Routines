# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 09:10:16 2013

@author: phil
"""

entlong=2**600000
ch_entlong=str(entlong)
sous_chaine="12345"
index=ch_entlong.find(sous_chaine)

print "Voici la position de ",sous_chaine," : ",index
