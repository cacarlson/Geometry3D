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
	print(compute_cut(cuts))
	return compute_cut(cuts)


_,_,_,_,cuts = construct_simplex(2/5)

alpha = 1/3
V,T,S,C, trims = construct_simplex(alpha)
# Jafar: Use trims for a candidate cut for testing purposes.
print(compute_cut([trim.plane for trim in trims]))

init = [x for cut in cuts for x in cut.plane.general_form()]
res = minimize(KCut, init)

print(res)

r = Renderer(backend='matplotlib')
r.add((T,'r',1),normal_length=0)
for s in S:
	r.add((s,'b',2),normal_length=0)
r.add((C,'g',3),normal_length=0)

data = res.x
opt_cut_data = [*zip(data[::4], data[1::4], data[2::4], data[3::4])]
opt_cuts = [Plane(*val) for val in opt_cut_data]
print(compute_cut(opt_cuts))
for cut in opt_cuts:
	r.add((intersection(cut, T), 'black',5),normal_length=0)
r.show()
