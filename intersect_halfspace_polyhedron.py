import math
import copy
import itertools
from Geometry3D import *

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
	pol = intersection(a,b)
	# pol is a convex polygon which has attributes such as its vertices and the plane it belongs to
	if not isinstance(pol, ConvexPolygon):		# check if pol is a polygon or not
			print("No cut ocurred!")		# we should be careful with this case
			return b
	s = Segment(pol.points[0], pol.points[1])					# first segment (edge) of pol
	c = pol.plane.n.cross(s.line.dv)							# vector parallel to a side of the cube on its certain face
	# Here pol.plane is the plane that pol belongs to and pol.plane.n is its unit normal vector
	length = 10 * max([seg.length() for seg in b.segment_set])	# length to ensure that the uncut part of b is captured ???
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
	return intersection(box, b)		# return the intersection


def compute_cut(cuts):
	# T is the simplex that we want to cut.
	# V is the list of vertices of T so that V[i] = vertex i of T.
	# S is the list of top simplices such that S[i] is the top simplex related to V[i] for i in {1,2,3,4}.
	# C is region of T \ union_of_(S[1], S[2], S[3], S[4]).
	# Here I assume these objects have already been defined globally.
	# cuts is a list of hyperplanes ordered so that cuts[i] cuts V[i].

	T_active = copy.deepcopy(T)		# This is the copy of the simplex that we will work with.
	cut_cost = float()				# cut cost is initially zero.
	for i,c in enumerate(cuts):
		if not isinstance(T_active, ConvexPolyhedron):
			raise TypeError("We should consider the case when T_active is not a polyhedron.")
		pol = intersection(c,T_active)		# pol is a polygon that cuts the "uncut" portion of the simplex T
		if not isinstance(pol,ConvexPolygon):
			continue			# we should be careful with this case
		cut_cost += grid_edges_cost(pol)		# compute the contribution of the grid edges (those parallel to the sides of T) to the cost of a cut.
		cut_cost += corner_edges_cost(pol)		# compute the contribution of the corner edges to the cost of a cut
		T_active = intersection(c,T_active,V[i])
	return cut_cost

def grid_edges_cost(a):
	'''
	Compute the contribution of the grid edges (those parallel to the sides of T) to the cost of a cut.
	'''
	# a is guaranteed to be a convex polygon
	grid_cost = float()
	a_grid = intersection(C, a)		# portion of a which contributes to grid edges' cost
	# check whether a_grid is a convex polygon.
	for seg in T.segment_set:
		grid_cost += a_grid.area() * seg.length() * abs(a_grid.plane.n * seg.line.dv.normalized())	# compute the contribution of edges parallel to seg, a side of T.
	return grid_cost

def corner_edges_cost(a):
	'''
	Compute the contribution of the corner edges to the cost of a cut
	'''
	# a is guaranteed to be a convex polygon
	corner_cost = float()
	for i,s in enumerate(S):
		a_corner = intersection(a, s)		# portion of a that contributes to corner edges' cost related to s
		if not isinstance(a_corner, ConvexPolygon):
			continue		# we should be careful with this case
		face_i = ConvexPolygon(tuple(s.point_set - {V[i]}))			# face of S[i] opposite to the vertex V[i]
		proj_vertices = set()										# This is the set of vertices of a projection of a_corner on face_i
		for p in a_corner.points:
			# We need to make sure that p is not equal to V[i] otherwise Line() function below will raise a value error.
			proj_vertex = intersection(Line(V[i], p),face_i)	# This is the "projection" of p on face_i
			if not isinstance(proj_vertex, Point):
				raise TypeError("Intersection is not a point.")				# We need to check whether p is of type Point.
			proj_vertices.add(proj_vertex)
		if len(proj_vertices) < 3:
			raise ValueError("To build a polygon the number of points cannot be less than 3")		# The case when projection of a_corner on face_i is not a polygon.
		proj_cut = ConvexPolygon(tuple(proj_vertices))
		corner_cost += proj_cut.area()				# compute the contribution of corner edges incident to V[i]. We should also include a multiplicative constant.
	return corner_cost

# The above libraries contain those that were used in the source code of the calc.intersection module
def construct_simplex(a):
	# The above libraries contain those that were used in the source code of the calc.intersection module

	# Setup the instance
	v1, v2, v3, v4 = Point(1,1,1), Point(-1,-1,1), Point(1,-1,-1), Point(-1,1,-1)	# Vertices of the simplex T ??
	V = [v1, v2, v3, v4]														# V is the list of vertices of T so that V[i] = vi of T.
	face1 = ConvexPolygon([v1, v2, v3])									# Faces of the T
	face2 = ConvexPolygon([v1, v2, v4])
	face3 = ConvexPolygon([v1, v3, v4])
	face4 = ConvexPolygon([v2, v3, v4])
	T = ConvexPolyhedron([face1, face2, face3, face4])						# Simplex that we want to cut

	# We need to define S = [s1, s2, s3, s4]
	S = []
	trims = []
	for v in V:
		pt_set = [v]
		for u in V:
			if u == v:
				continue
			pt_set.append(Point((1-a)*v.pv() + (a)*u.pv()))
		faces = []
		for face in itertools.combinations(pt_set, 3):
			faces.append(ConvexPolygon(face))
		S.append(ConvexPolyhedron(faces))
		trims.append(ConvexPolygon(list(pt_set)[1:]))

	# We need to define C= T \ union_of_(S[1], S[2], S[3], S[4])
	C = copy.deepcopy(T)
	# for i, trim in enumerate(trims):
	# 	print(V[i],trim.plane)
	# 	C = (inter_halfspace_convexpolyhedron(trim.plane, C, V[i]))

	return (V, T, S, C)
