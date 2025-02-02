import math
import copy
import itertools
from Geometry3D import *
import matplotlib.pyplot as plt

side_length =  2 * (2 ** 0.5) ## the simplex side length is sqrt(8) in the current implementation
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

	# pol could be nothing, point or line (if it isn't a polygon)
	# if this is the case, then either the intersection is pol or it is b
	if not isinstance(pol, ConvexPolygon):		# check if pol is a polygon or not
		cp = b._get_center_point()
		if ((cp.pv() - a.p.pv())*a.n) > 0:
			return b
		return pol
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
	u1, u2, u3, u4 = Point(v1.pv() + e), Point(v2.pv() + e), Point(v3.pv() + e), Point(v4.pv() + e)	# second group of four vertices on the opposite face of the cube
	f1 = ConvexPolygon([v1, v2, v3, v4])					# faces 1 through 6 of the cube
	f2 = ConvexPolygon([u1, u2, u3, u4])
	f3 = ConvexPolygon([u1, u2, v1, v2])
	f4 = ConvexPolygon([v3, v4, u3, u4])
	f5 = ConvexPolygon([u1, u3, v1, v3])
	f6 = ConvexPolygon([v2, v4, u2, u4])

	box = ConvexPolyhedron([f1, f2, f3, f4, f5, f6])		# cube that mimics the halfspace

	# r = Renderer(backend='matplotlib')
	# r.add((b,'r',1),normal_length=0)
	# r.add((box,'b',1),normal_length=0)
	# r.show()
	return intersection(box, b)		# return the intersection

def compute_cut(cuts, T, V, S, C, alpha):
	# T is the simplex that we want to cut.
	# V is the list of vertices of T so that V[i] = vertex i of T.
	# S is the list of top simplices such that S[i] is the top simplex related to V[i] for i in {1,2,3,4}.
	# C is region of T \ union_of_(S[1], S[2], S[3], S[4]).
	# cuts is a list of hyperplanes ordered so that cuts[i] cuts V[i].

	faces = T.convex_polygons #[(ConvexPolygon(tuple(s.point_set - {V[i]})),i) for i,s in enumerate(S)]
	#faces_active = [(ConvexPolygon(tuple(s.point_set - {V[i]})),i) for i,s in enumerate(S)]
	#cut_cost_2d = float()				# cut cost is initially zero.
	grid_cost_2d = float()
	corner_cost_2d = float()

	T_active = copy.deepcopy(T)		# This is the copy of the simplex that we will work with.
	#cut_cost_3d = float()				# cut cost is initially zero.
	grid_cost_3d = float()
	corner_cost_3d = float()

	min_dis = next(iter(T.segment_set)).length()

	for i,c in enumerate(cuts):
		#print("Cut: \n")

		for v in V:
			min_dis = min(min_dis, distance(v,c))

		if not isinstance(T_active, ConvexPolyhedron):
			raise TypeError("We should consider the case when T_active is not a polyhedron.")

		# r = Renderer(backend='matplotlib')
		# #r.add((c,'r',1),normal_length=0)
		# r.add((T_active,'b',1),normal_length=0)
		# r.show()

		pol = intersection(c,T_active)		# pol is a polygon that cuts the "uncut" portion of the simplex T
		# if (isinstance(pol,ConvexPolygon)):
		# 	print("Area: ", pol.area())
		if not isinstance(pol,ConvexPolygon):
			continue

		# # lets compute max point distance in pol
		# mx_seg = 0
		#
		# for segs in pol.segments():
		# 	mx_seg = max(mx_seg, segs.length())
		#
		# if mx_seg < .001:
		# 	continue	# we should be careful with this case
		# 				# Note that we aren't adjusting T_active but the only way this should
		# 				# come up is if we are removing a small volume which shouldn't be the case

		grid_cost_3d += grid_edges_cost_3d(pol, T, C)		# compute the contribution of the grid edges (those parallel to the sides of T) to the cost of a cut.
		corner_cost_3d += corner_edges_cost_3d(pol, S, V)		# compute the contribution of the corner edges to the cost of a cut

		# each face of T
		for face in faces:
			#print("Cut: ", c)

			#print("Face: ", face)
			face_active = copy.deepcopy(face)
			face_active = intersection(face_active, T_active)
			#print("Face Active: ", face_active)

			S_faces = [intersection(face, s) for s in S]
			S_faces = list(filter(lambda s: s != None, S_faces))
			C_face = intersection(C, face)
			verts_face = face.points

			if not isinstance(face_active, ConvexPolygon):
				#print("Face not active!")
				# if None then pruned face
				# do we care about segments or points?
				# raise TypeError("We should consider the case when T_active is not a polyhedron.")
				continue

			#local_cut = intersection(c, T)

			face_pol = intersection(c, face_active)
			#print("Face pol: ", face_pol)

			#r = Renderer(backend='matplotlib')
			#
			# r.add((T,'r',1),normal_length=0)
			# # for s in S:
			# # 	r.add((s,'b',2),normal_length=0)
			# # r.add((C,'g',3),normal_length=0)
			# #r.add((local_cut, 'green', 5), normal_length=0)
			# #r.add((intersection(c, T), 'green', 5), normal_length=0)
			# r.add((face_active, 'blue', 5), normal_length=0)
			# r.add((face_pol, 'black', 5), normal_length=0)
			# # for cut in cuts:
			# # 	r.add((intersection(cut, T), 'black',5),normal_length=0)
			# r.show()

			if not isinstance(face_pol, Segment):
				#print("Intersect not segment")
				# either face was removed
				# or cut cuts along same as previous cut and no cost
				continue
			grid_cost_2d += grid_edges_cost_2d(face_pol, face, C_face)		# compute the contribution of the grid edges (those parallel to the sides of T) to the cost of a cut.
			#print("Cont: ", corner_edges_cost_2d(face_pol, S_faces, verts_face))
			# return 0
			corner_cost_2d += corner_edges_cost_2d(face_pol, S_faces, verts_face)

		T_active = inter_halfspace_convexpolyhedron(c,T_active,V[i])

#	print("Min Dis: ", min_dis)
	pen = 2**(6*(2/3-min(min_dis, 2*alpha)))

	for seg in T.segment_set:
		if seg in T_active:
			pen = max(pen,2**(3))
			break

	# for v in V:
	# 	min_d = 1/3
	# 	for point in T_active.point_set:
	# 		#print("distance: ", distance(point,v)/2)
	# 		min_d = min(min_d, distance(point,v)/2)
	# 	pen = max(pen, 2**(9*(1/3-min_d)))
	# print("min distance: ", min_d)
	# print("pen: ", pen)
#	print("3d Cost: ", corner_cost_3d +grid_cost_3d)
#	print("2d Cost: ", corner_cost_2d +grid_cost_2d)

##	pen = 1 # for debugging
	return pen*(corner_cost_2d +grid_cost_2d), pen*(grid_cost_3d + corner_cost_3d) 

def grid_edges_cost_3d(a, T, C):
	'''
	Compute the contribution of the grid edges (those parallel to the sides of T) to the cost of a cut.
	'''
	# a is guaranteed to be a convex polygon
	grid_cost = float()

	a_grid = intersection(C, a)		# portion of a which contributes to grid edges' cost

	if (a_grid == None):
		#print ("Grid Cost: ", grid_cost)
		return grid_cost
	# if not isinstance(a_grid, ConvexPolyhedron):
	# 	return grid_cost
	# check whether a_grid is a convex polygon.
	for seg in T.segment_set:
		# we don't have to multiply by side length (be consistent)
		# LP cost depends on this choice
		# Jafar: Removed seg.length() from the below product
		# note that in this implementation seg.length = sqrt(8)
		#the formula below should work even if we rescale the simplex
		grid_cost += (2 * (2 ** 0.5) / (side_length ** 2)) * a_grid.area() * abs(a_grid.plane.n * seg.line.dv.normalized())	# compute the contribution of edges parallel to seg, a side of T.
	return grid_cost

def corner_edges_cost_3d(a, S, V):
	'''
	Compute the contribution of the corner edges to the cost of a cut
	'''
	# a is guaranteed to be a convex polygon
	corner_cost = float()
	penalty = 0
	for i,s in enumerate(S):
		#print(a)
		#print(s)
		a_corner = intersection(a, s)		# portion of a that contributes to corner edges' cost related to s
		# TODO: What if None?
		if not isinstance(a_corner, ConvexPolygon):
			# this could be empty
			continue		# we should be careful with this case
		# TODO: Maybe we keep this around?
		face_i = ConvexPolygon(tuple(s.point_set - {V[i]}))			# face of S[i] opposite to the vertex V[i]
		proj_vertices = set()
										# This is the set of vertices of a projection of a_corner on face_i
		for p in a_corner.points:
			# We need to make sure that p is not equal to V[i] otherwise Line() function below will raise a value error.
			# dist = distance(V[i], p)
			# if dist < 1/30:
			# 	penalty += 100 * (1/30 - dist)
			proj_vertex = intersection(Line(V[i], p),face_i)	# This is the "projection" of p on face_i
			if not isinstance(proj_vertex, Point):
				# Note to Charlie: This should not happen so error works
				raise TypeError("Intersection is not a point.")				# We need to check whether p is of type Point.
			proj_vertices.add(proj_vertex)
		if len(proj_vertices) < 3:
			# Note to Charlie: This should not happen either
			raise ValueError("To build a polygon the number of points cannot be less than 3")		# The case when projection of a_corner on face_i is not a polygon.
		proj_cut = ConvexPolygon(tuple(proj_vertices))
		# This seems off:
		#corner_cost += (2)*(proj_cut.area() * math.sqrt(3) / 2)	# Jafar: Included the multiplicative constant. Removed side length contribution from both types of cuts.
		corner_cost += (4 * (3 ** 0.5) / side_length ** 2 ) * proj_cut.area() # + penalty

		#*edge weights (area of proj_cut) / (area of face_i ) * (Cost of face_i)
		# cost of face_i should be the same as moving it a little above
		# area of triangle times sum of unit vectors times unit vectors
		# root(6)/2
		# root(3)/2 * side length <-
		# compute the contribution of corner edges incident to V[i]. We should also include a multiplicative constant.
	#print("Corner Cost: ", corner_cost)

	return corner_cost

def grid_edges_cost_2d(a, face, c_face):
	'''
	Compute the contribution of the grid edges (those parallel to the sides of T) to the cost of a cut.
	'''
	# a is guaranteed to be a convex polygon
	grid_cost = float()

	a_grid = intersection(c_face, a)	# portion of a which contributes to grid edges' cost

	if not isinstance(a_grid, Segment):
		return grid_cost
    
	segs = [seg.line.dv.normalized() for seg in face.segments()]
	altitudes = [segs[0] + segs[1] if segs[0] * segs[1] > 0 else segs[0] - segs[1],
	             segs[0] + segs[2] if segs[0] * segs[2] > 0 else segs[0] - segs[2],
				 segs[1] + segs[2] if segs[1] * segs[2] > 0 else segs[1] - segs[2]
	]

	for h in altitudes:
		grid_cost = max(grid_cost, (4/(3 ** 0.5)) * a_grid.length() / side_length * abs(a_grid.line.dv.normalized() * h.normalized()))

	return grid_cost

def corner_edges_cost_2d(a, S_face, verts_face):
	'''
	Compute the contribution of the corner edges to the cost of a cut
	'''

	# print(a)
	# print(S_face)
	# print(verts_face)
	# print("\n")

	# a is guaranteed to be a convex polygon
	corner_cost = 0

	for s in S_face:
		a_corner = intersection(a, s)		# portion of a that contributes to corner edges' cost related to s
		# TODO: What if None?
		if not isinstance(a_corner, Segment):
			# this could be empty
			continue		# we should be careful with this case
		# TODO: Maybe we keep this around?

		# Find original vert in face
		for v in verts_face:
			if v in s:
				#print("found")
				orgin = v
				seg = Segment(*tuple(set(s.points) - {v}))

		#seg_i = Segment(*tuple(set(s.points) - {verts[i]}))			# face of S[i] opposite to the vertex V[i]
		proj_vertices = []
										# This is the set of vertices of a projection of a_corner on face_i
		for p in [a_corner.start_point, a_corner.end_point]:
			# We need to make sure that p is not equal to V[i] otherwise Line() function below will raise a value error.
			# dist = distance(V[i], p)
			# if dist < 1/30:
			# 	penalty += 100 * (1/30 - dist)
			proj_vertex = intersection(Line(orgin, p),seg)	# This is the "projection" of p on face_i
			if not isinstance(proj_vertex, Point):
				# Note to Charlie: This should not happen so error works
				raise TypeError("Intersection is not a point.")				# We need to check whether p is of type Point.
			proj_vertices.append(proj_vertex)
		if len(proj_vertices) != 2:
			# Note to Charlie: This should not happen either
			raise ValueError("To build a polygon the number of points cannot be less than 3")		# The case when projection of a_corner on face_i is not a polygon.
		proj_cut = Segment(proj_vertices[0], proj_vertices[1])
		# This seems off:
		#corner_cost += (2)*(proj_cut.area() * math.sqrt(3) / 2)	# Jafar: Included the multiplicative constant. Removed side length contribution from both types of cuts.
		corner_cost += 2 * proj_cut.length() / side_length # + penalty

		#*edge weights (area of proj_cut) / (area of face_i ) * (Cost of face_i)
		# cost of face_i should be the same as moving it a little above
		# area of triangle times sum of unit vectors times unit vectors
		# root(6)/2
		# root(3)/2 * side length <-
		# compute the contribution of corner edges incident to V[i]. We should also include a multiplicative constant.
	#print("Corner Cost: ", corner_cost)

	return corner_cost

# Add 1D cost function for each edge
def construct_simplex(a):
	# The above libraries contain those that were used in the source code of the calc.intersection module
	#global V, T, C, S, trims	# Jafar: These variables are accessed by compute_cut(), grid_edges_cost(), corner_edges_cost() functions.
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
			faces.append(ConvexPolygon(tuple(face)))			# Jafar: Input to ConvexPolygon() should be a tuple
		S.append(ConvexPolyhedron(tuple(faces)))				# Jafar: Input to ConvexPolyhedron() should be a tuple
		trims.append(ConvexPolygon(tuple(pt_set[1:])))			# Jafar: Input to ConvexPolygon() should be a tuple

	# We need to define C= T - union_of_(S[1], S[2], S[3], S[4])
	C = copy.deepcopy(T)
	for i, trim in enumerate(trims):
		#print(V[i],trim.plane)
		C = (inter_halfspace_convexpolyhedron(trim.plane, C, V[i]))

	return (V, T, S, C, trims)		#Jafar: Added trims to the output for testing purposes.

def LP_cost(alpha):
	'''
	Compute the LP cost for the given value of alpha. The function returns the 2d and 3d costs.
	'''	
	##LP cost: (d+1)/2 * (1 + (d - 1) * alpha ** d)
	## we have 4 copies of the 2d instance
	return (6 * (1 + alpha ** 2), 2 * (1 + 2 * alpha ** 3))