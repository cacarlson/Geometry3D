import math
import copy
from Geometry3D import *

import numpy as np
import scipy as sp
from scipy.optimize import minimize

from intersect_halfspace_polyhedron import *

# Globals .... okay for now
def KCut(data):
	cut_data = [*zip(data[::4], data[1::4], data[2::4], data[3::4])]
	#print(cut_data)
	cuts = [Plane(*val) for val in cut_data]
	#print(cuts)
	return compute_cut(cuts)

alpha = 1/3
V,T,S,C, trims = construct_simplex(alpha)		# Jafar: Use trims for a candidate cut for testing purposes.
print(compute_cut([trim.plane for trim in trims]))

init = [x for trim in trims for x in trim.plane.general_form()]
data = init
cut_data = [*zip(data[::4], data[1::4], data[2::4], data[3::4])]
cuts = [Plane(*val) for val in cut_data]
print(compute_cut(cuts))


res = minimize(KCut, init)
print(res)

r = Renderer(backend='matplotlib')
r.add((T,'r',1),normal_length=0)
# for s in S:
# 	r.add((s,'b',2),normal_length=0)
# r.add((C,'g',3),normal_length=0)

data = res.x
cut_data = [*zip(data[::4], data[1::4], data[2::4], data[3::4])]
cuts = [Plane(*val) for val in cut_data]
print(compute_cut(cuts))
for cut in cuts:
	r.add((intersection(cut, T), 'black',5),normal_length=0)
r.show()
