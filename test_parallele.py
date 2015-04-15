# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 13:33:59 2012

@author: dbonnet
"""

#import lib_parallele
#
#ABC = lib_parallele.obj(2,3)

from IPython import parallel as p
c = p.Client()
dview = c[:]
e0 = c[0]
dview.block = True

#chemin de la librairie a importer
PATH = '/home/dbonnet/Bureau'

with dview.sync_imports():
    """
    import de sys sur tous les noeuds
    """
    import sys
    
def append_sys_path(path):
    """
    pour l'appliquer sur tous les noeuds
    """
    sys.path.append(path)
    
dview.apply(append_sys_path, PATH)

with dview.sync_imports():
    """
    Maintenant que le sys.path est a jour sur tous les noeuds, on importe la
    librairie
    """
    import lib_parallele
    import numpy


def Test_passage_numpy(B):
    """
    Un test pour envoyer ou recevoir un tableau numpy
    (impossible de passer un tableau vide)
    """
    return B
    
def Test(A):
    """
    Test pour passer un objet
    """
    return A.plus()    

def Truc(a, b):
    """
    test pour crer un objet a partir de la librairie importee
    """
    return lib_parallele.obj(a, b).plus(), numpy.array([a, b])
    
def func():
    """
    test d'une fonction
    """
    return numpy.array([2])

def fun():
    """
    test d'une fonction
    """
    return 2+2

#map reparti l'operation Truc(a, b) pour 10 couples d'arguments sur les noeuds
res2 = dview.map(Truc, range(10), range(10))

#on cree un objet qu'on passe
a = lib_parallele.obj(3, 4)
res3 = dview.map(Test, [a]*10)

#apply realise la fonction Test avec l'argument a sur chaque coeur
res4 = dview.apply(Test, a)

res5 = dview.apply(func)

res6 = dview.apply(fun)

res7 = dview.apply(Test_passage_numpy, [])

# On cree un objet sur chaque coeur
res8 = dview.apply(lib_parallele.obj, 3, 2)

#print res
print 'res2 : ', res2
print 'res3 : ', res3
print 'res4 : ', res4
print 'res5 : ', res5
print 'res6 : ', res6
print 'res7 : ', res7
print 'res8[0].a :', res8[0].a