import math
import copy
from Geometry3D import *

import numpy as np
import scipy as sp
from scipy.optimize import minimize

from intersect_halfspace_polyhedron import *

import matplotlib.pyplot as plt

# Globals .... okay for now
def KCut(data):
	cut_data = [*zip(data[::4], data[1::4], data[2::4], data[3::4])]
	cuts = [Plane(*val) for val in cut_data]
	v = compute_cut(cuts, T, V, S, C)[0]
	print(v)
	return v


set_eps(1e-10)

# Jafar: Use trims for a candidate cut for testing purposes.
# print(compute_cut([trim.plane for trim in trims]))
# r = Renderer(backend='matplotlib')
# r.add((T,'r',1),normal_length=0)
# for s in S:
# 	r.add((s,'b',2),normal_length=0)
# r.add((C,'g',3),normal_length=0)
#
# for cut in cuts:
# 	r.add((intersection(cut, T), 'black',5),normal_length=0)
# r.show()

#init = [x for cut in cuts for x in cut.plane.general_form()]

init = [3,-2,4,4,
		-1,-2,5,3.5,
		2,1,-3,3,
		-4,0,-5,2]
# init = [ 2.8678988,  -2.0673227,   3.67714101,  4.38827445, -0.87509477, -1.72442337,
#   			4.37484964,  4.5028203,   1.86045445,  1.07570424, -2.64171212,  3.42644829,
#  			-3.89068623,  0.41426979, -4.7023626,   2.85373303]
# init = [ 2.90267946, -2.18060441,  3.59752785,  4.31949062, -0.74468934, -1.53388329,
#  			3.86423413,  5.26364051,  1.93588653,  1.15961594, -2.67864448,  3.21193694,
#  			-3.60791907,  0.84149244, -4.10121173,  4.21413749]
# init = [ 2.86789808, -2.06732361,  3.67714059,  4.388274,   -0.8750945,  -1.72442267,
#   	4.37484828,  4.50282188,  1.86045413,  1.07570371, -2.64171244,  3.42644778,
#  	-3.8906908,   0.41427025, -4.70236062,  2.85373489]

#rint(compute_cut([trim.plane for trim in trims]))
#init_data = [*zip(init[::4], init[1::4], init[2::4], init[3::4])]
#init_cuts = [Plane(*val) for val in init_data]

alpha = 1/3
# cut_costs = []
# corner_costs = []
# grid_costs = []
# points = [alpha + off/6 for off in np.linspace(-.6,6/2-6/3-.01,20)]
# for beta in points:

_,_,_,_,cuts = construct_simplex(0.2)
V,T,S,C,trims = construct_simplex(alpha)

#init = [x for cut in cuts for x in cut.plane.general_form()]
init_data = [*zip(init[::4], init[1::4], init[2::4], init[3::4])]
init_cuts = [Plane(*val) for val in init_data]

# 	print(beta)
# 	cut_cost, grid_cost, corner_cost = compute_cut(init_cuts)
# 	cut_costs.append(cut_cost) #/beta**2/3)
# 	grid_costs.append(grid_cost)
# 	corner_costs.append(corner_cost)

# 	#print(compute_cut(init_cuts)/beta**2/3)
# plt.plot(points, cut_costs, 'r') # plotting t, a separately
# plt.plot(points, grid_costs, 'b') # plotting t, b separately
# plt.plot(points, corner_costs, 'g') # plotting t, c separately
# plt.show()

# r = Renderer(backend='matplotlib')
# r.add((T,'r',1),normal_length=0)
# for s in S:
# 	r.add((s,'b',2),normal_length=0)
# r.add((C,'g',3),normal_length=0)

# for cut in init_cuts:
# 	print(cut)
# 	r.add((intersection(cut, T), 'black',5),normal_length=0)
# r.show()


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
print(compute_cut(opt_cuts, T, V, S, C))
for cut in opt_cuts:
 	r.add((intersection(cut, T), 'black',5),normal_length=0)
r.show()
