import math
import copy
from Geometry3D import *
# The above libraries contain those that were used in the source code of the calc.intersection module

def inter_halfspace_convexpolyhedron(a,b,v):
	'''
	I used the idea that we discussed couple of meetings ago: defining a large enough cube so that it contains
	the "uncut" polyhedron by using the intersection-of-two-polyhedrons function from the library.
	'''
	# a is a plane
	# b is a convex polyhedron
	# v is the vertex that we want to cut
	a = copy.deepcopy(a)
	b = copy.deepcopy(b)
	v = copy.deepcopy(v)

	if not isinstance(a, Plane) or not isinstance(b, ConvexPolyhedron) or not isinstance(v, Point):
		raise TypeError("Incorrect Input")
	pol = inter_plane_convexpolyhedron(a,b)
	# pol is a convex polygon which has attributes such as its vertices and the plane it belongs to
	if not isinstance(pol, ConvexPolygon):		# check if pol is a polygon or not
			print("No cut ocurred!")		# we should be careful with this case
			return b
	s = Segment(pol.points[0], pol.points[1])					# first segment (edge) of pol
	c = pol.plane.n.cross(s.line.dv)							# vector parallel to a side of the cube on its certain face
	# Here pol.plane is the plane that pol belongs to and pol.plane.n is its unit normal vector
	length = 10 * max([seg.length() for seg in b.segmen_set])	# length to ensure that the uncut part of b is captured ???
	# We want to build a cube of side-length equal to 2*length
	c = c.normalized() * length
	d = s.line.dv.normalized() * (length - s.length() / 2)		# vector parallel to the other side of the cube on the same face
	x, y = s.start_point.pv() - d, s.end_point.pv() + d
	v1, v2, v3, v4 = Point(x + c), Point(x - c), Point(y + c), Point(y - c)	# first group of four vertices on the face of the cube
	e = pol.plane.n * 2 * length
	if pol.plane.n * v.pv() > pol.plane.n * pol.plane.p.pv():	# we decide for the "direction" of the halfspace (cube)
			e = e * (-1)
	u1, u2, u3, u4 = Point(v1 + e), Point(v2 + e), Point(v3 + e), Point(v4 + e)	# second group of four vertices on the opposite face of the cube
	f1 = ConvexPolygon(tuple(v1, v2, v3, v4))					# faces 1 through 6 of the cube
	f2 = ConvexPolygon(tuple(u1, u2, u3, u4))
	f3 = ConvexPolygon(tuple(u1, u2, v1, v2))
	f4 = ConvexPolygon(tuple(v3, v4, u3, u4))
	f5 = ConvexPolygon(tuple(u1, u3, v1, v3))
	f6 = ConvexPolygon(tuple(v2, v4, u2, u4))
	box = ConvexPolyhedron(tuple(f1, f2, f3, f4, f5, f6))		# cube that mimics the halfspace
	return inter_convexpolyhedron_convexpolyhedron(box, b)		# return the intersection


def compute_cut(cuts):
	# T is the simplex that we want to cut.
	# Si is the top simplex related to ei, that is, vertex i of T for i in {1,2,3,4}.
	# C is region of T \ union_of_(S1, S2, S3, S4).
	# Here I assume these objects have already been defined globally.
	# cuts is a list of hyperplanes ordered so that cuts[i] cuts ei.

	T_active = copy.deepcopy(T)		# This is the copy of the simplex that we will work with.
	cut_cost = float()				# cut cost is initially zero.
	for i,c in enumerate(cuts):
		if not isinstance(T_active, ConvexPolyhedron):
			raise TypeError("We should consider the case when T_active is not a polyhedron.")
		pol = inter_plane_convexpolyhedron(c,T_active)		# pol is a polygon that cuts the "uncut" portion of the simplex T
		if not isinstance(pol,ConvexPolygon):
			continue			# we should be careful with this case
		cut_cost += grid_edges_cost(pol)		# compute the contribution of the grid edges (those parallel to the sides of T) to the cost of a cut.
		cut_cost += corner_edges_cost(pol)		# compute the contribution of the corner edges to the cost of a cut
		T_active = inter_halfspace_convexpolyhedron(c,T_active,E[i])	# E is the list of vertices of T so that E[i] = ei.
	return cut_cost

def grid_edges_cost(a):
	'''
	Compute the contribution of the grid edges (those parallel to the sides of T) to the cost of a cut.
	'''
	# a is guaranteed to be a convex polygon
	grid_cost = float()
	a_grid = inter_convexpolygon_convexPolyhedron(C, a)		# portion of a which contributes to grid edges' cost
	# check whether a_grid is a convex polygon.
	for seg in T.segment_set:
		grid_cost += a_grid.area() * seg.length() * a_grid.plane.n * seg.line.dv.normalized()	# compute the contribution of edges parallel to seg, a side of T.
	return grid_cost

def corner_edges_cost(a):
	'''
	Compute the contribution of the corner edges to the cost of a cut
	'''
	# a is guaranteed to be a convex polygon
	corner_cost = float()
	for x in {S1,S2,S3,S4}:
		a_corner = inter_convexpolygon_convexPolyhedron(a, x)		# portion of a which contributes to corner edges' cost
		#........
		#Not sure how to compute the contribution of corner edges to the cut
	return corner_cost
