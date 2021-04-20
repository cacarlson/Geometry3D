import math
import copy
from Geometry3D import *

import numpy as np
import scipy as sp
from scipy.optimize import minimize

from intersect_halfspace_polyhedron import *

# Globals .... okay for now
def KCutSimplex(args):
	print(args)
	x,y = args
	return abs(x+y)

alpha = 1/3
V,T,S,C, trims = construct_simplex(alpha)		# Jafar: Use trims for a candidate cut for testing purposes.
print(compute_cut([trim.plane for trim in trims]), alpha)
r = Renderer(backend='matplotlib')
r.add((T,'r',1),normal_length=0)
for s in S:
	r.add((s,'b',2),normal_length=0)
r.add((C,'g',3),normal_length=0)
r.show()
