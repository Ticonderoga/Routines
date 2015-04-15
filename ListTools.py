# -*- coding: iso-8859-15 -*-
def union(L1,L2) :
	"""list1=range(5) # [0,1,2,3,4]
	list2=range(3,7) # [3,4,5,6]
	union(list1,list2) = [0,1,2,3,4,5,6]"""
	return L1+filter(lambda x:x not in L1,L2)


def intersect(L1,L2) :
	"""list1=range(5) # [0,1,2,3,4]
	list2=range(3,7) # [3,4,5,6]
	intersection(list1,list2)=[3,4]"""
	return filter(lambda x:x in L1,L2)

def difference(L1,L2) :
	"""list1=range(5) # [0,1,2,3,4]
	list2=range(3,7) # [3,4,5,6]
	difference(list1,list2)=[0,1,2]"""
	return filter(lambda x:x not in L2,L1)

def distinct(L1,L2) :
	"""list1=range(5) # [0,1,2,3,4]
	list2=range(3,7) # [3,4,5,6]
	distinct(list1,list2)=[0,1,2,5,6]"""
	return filter(lambda x:x not in L2,L1)+filter(lambda x:x not in L1,L2)

if __name__=='__main__' :
	list1=range(5)
	list2=range(3,7)
	print "list1 :",list1
	print "list2 :",list2
	print "union l1,l2 :",union(list1,list2)
	print "intersect l1,l2 :",intersect(list1,list2)
	print "difference l1,l2 :",difference(list1,list2)
	print "distinct l1,l2 :",distinct(list1,list2)
	
