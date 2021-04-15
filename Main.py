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

V,T,S,C = construct_simplex(1/3)
r = Renderer(backend='matplotlib')
r.add((T,'r',1),normal_length=0)
for s in S:
	r.add((s,'b',2),normal_length=0)
r.add((C,'g',3),normal_length=0)
r.show()
