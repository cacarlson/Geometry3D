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
    cuts = [Plane(*val) for val in cut_data]
    print(data)
    computed_cuts = compute_cut(cuts, T, V, S, C)
    print(computed_cuts)
    return computed_cuts


set_eps(1e-10)

## alpha = 1/3
for alpha in [1/3 - 0.1, 1/3 -0.05, 1/3 - 0.01, 1/3 + 0.001, 1/3 + 0.05]: ##[1/3 + x / 300.0 for x in range(1, 10)]:
    _,_,_,_,cuts  = construct_simplex(alpha)
    V,T,S,C, trims = construct_simplex(1/3)
    init = [x for cut in cuts for x in cut.plane.general_form()]
    init_data = [*zip(init[::4], init[1::4], init[2::4], init[3::4])]
    init_cuts = [Plane(*val) for val in init_data]

    cost,_,_ = compute_cut(init_cuts, T, V, S, C)
    print(alpha - 1/3, ":   ", cost, cost / (9 * alpha ** 2))

# -0.1 :    1.333333333333334 2.7210884353741513
# -0.04999999999999999 :    1.3333333333333344 1.8454440599769335
# -0.010000000000000009 :    1.333333333333334 1.4170829347787588

# 0.0010000000000000009 :    1.2646325081505743 1.2570787221094188
# 0.04999999999999999 :    1.6624866099897035 1.2570787221094168

    # r = Renderer(backend='matplotlib')
    # r.add((T,'r',1),normal_length=0)
    # for s in S:
    #     r.add((s,'b',2),normal_length=0)
    # r.add((C,'g',3),normal_length=0)
    # r.show()

exit()

for cut in init_cuts:
    print(cut)
    r.add((intersection(cut, T), 'black',5),normal_length=0)
r.show()
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
