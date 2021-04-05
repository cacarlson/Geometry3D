from Geometry3D import *
import copy

pa = Point(1,1,1)
pb = Point(-1,-1,1)
pc = Point(-1,1,-1)
pd = Point(1,-1,-1)

f1 = ConvexPolygon((pa, pb, pc))
f2 = ConvexPolygon((pa, pb, pd))
f3 = ConvexPolygon((pa, pc, pd))
f4 = ConvexPolygon((pb, pc, pd))

cpg = ConvexPolyhedron((f1,f2,f3,f4))

#h = copy.deepcopy(f1.plane)
r = Renderer(backend='matplotlib')


h = HalfSpace(f1.plane)
h.move(Vector(1,-1,1))

#print("Inter of: ", cpg, h)
cpg = intersection(cpg, h)
#print(cpg)


r.add((cpg,'b',2),normal_length=0)
r.add((f1.move(Vector(1,-1,1)),'r',2),normal_length=0)
r.show()
