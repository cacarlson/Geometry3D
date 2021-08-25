import math
import copy
from Geometry3D import *

import numpy as np
import scipy as sp
from scipy.optimize import minimize

from intersect_halfspace_polyhedron import *

import matplotlib.pyplot as plt

# Globals .... okay for now
# T is the simplex that we want to cut.
T = None
# V is the list of vertices of T so that V[i] = vertex i of T.
V = None
# S is the list of top simplices such that S[i] is the top simplex related to V[i] for i in {1,2,3,4}.
S = None
# C is region of T \ union_of_(S[1], S[2], S[3], S[4]).
C = None
# trim is the planes of the top simplex regions
trims = None

def KCut(data, alpha, weight2d):
    cut_data = [*zip(data[::4], data[1::4], data[2::4], data[3::4])]    
    cuts = [Plane(*val) for val in cut_data]
    try: 
        costs = compute_cut(cuts, T, V, S, C, alpha)
        print("value: ", weight2d * costs[0] + costs[1])
        return (weight2d * costs[0] + costs[1])
    except:
        print("!!! Error/exception")
        return(10)


def OptimizeCut(hyper_params):
    alpha = hyper_params[0]
    weight2d = hyper_params[1]
    _,_,_,_,cuts = construct_simplex(alpha+0.05)
    init = [x for cut in cuts for x in cut.plane.general_form()]

    global V, T, C, S, trims
    V,T,S,C,trims = construct_simplex(alpha)    
    cut = minimize(KCut, init, args=(alpha, weight2d), tol=0.001, method='Powell')
    LP = LP_cost(alpha)

    print("=== DEBUG ===")
    print("integral: ", cut.fun)
    print("fractional: ", LP[0] * weight2d + LP[1])
    print("=============")
 
    
    return cut.fun / (LP[0] * weight2d + LP[1])

def main():   
    bestAlpha = 0.0
    bestWeight = 0.0
    bestGap = 0
    # alpha = 0.43333333333333335
    # weight2d = 0.21052631578947367
    # gap = OptimizeCut([alpha, weight2d])
    # print("Gap: ", gap, "for ", (alpha, weight2d))


    for alpha in np.linspace(0.34,0.48,10):
        for weight2d in np.linspace(0.0,4,20):   
            gap = OptimizeCut([alpha, weight2d])
            if gap > bestGap:
                print("Gap: ", gap, "for ", (alpha, weight2d))
                bestGap = gap
                bestWeight = weight2d
                bestAlpha = alpha

    print("Gap:", bestGap, "Best parameters:", (bestAlpha, bestWeight))

set_eps(1e-10)
if __name__ == "__main__":
    main()
