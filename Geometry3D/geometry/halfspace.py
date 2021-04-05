# -*- coding: utf-8 -*-
"""HalfSpace Module"""
from .body import GeoBody
from .point import Point
from .plane import Plane
from .polygon import ConvexPolygon
from .polyhedron import ConvexPolyhedron
from ..utils.vector import Vector
from .line import Line
from ..utils.constant import get_eps
from ..utils.logger import get_main_logger
from .segment import Segment
import math
import copy

class HalfSpace(GeoBody):
    """
    **Input:**

    - HalfLine(Plane, Vector)
    """
    class_level = 7 # the class level of HalfLine
    def __init__(self,pln):
        self.plane = copy.deepcopy(pln)

    def __eq__(self,other):
        return (self.plane == other.plane)

    def __repr__(self):
        return "HalfSpace({})".format(self.plane)

    def __contains__(self, other):
        """Checks if other lies in halfspace"""
        if isinstance(other,Point):
            return (((other.pv() - self.plane.p)*self.plane.n) <= 0)
        elif isinstance(other,Segment):
            return (other.start_point in self) and (other.end_point in self)
        elif isinstance(other,ConvexPolygon):
            for point in other.points:
                if not point in self:
                    return False
            return True
        elif isinstance(other,ConvexPolyhedron):
            for plg in other.convex_polygons:
                if not plg in self:
                    return False
            return True
        else:
            get_main_logger().warning("Calling type {} in type {} which is always False".format(type(other),type(self)))
            return False

    def __hash__(self):
        """return the hash of a HalfSpace"""
        return hash(("HalfSpace",round(self.n[0],SIG_FIGURES),round(self.n[1],SIG_FIGURES),round(self.n[2],SIG_FIGURES),round(self.n * self.p.pv(),SIG_FIGURES)))

    def move(self,v):
        """Return the HalfSpace that you get when you move self by vector v, self is also moved"""
        if isinstance(v,Vector):
            self.plane.move(v)
            return HalfSpace(self.plane)
        else:
            return NotImplementedError("The second parameter for move function must be Vector")

    def parametric(self):
        """Returns (u, v, w) so that you can build the equation
           _   _    _    _
        E: x <= u + rv + sw ; (r, s) e R

        to describe the plane (a point and two vectors).
        """
        s = solve([list(self.n) + [0]])
        # Pick a first vector orthogonal to the normal vector
        # there are infinitely many solutions, varying in direction
        # and length, so just choose some values
        v = Vector(*s(1, 1))
        assert v.orthogonal(self.n)
        # Pick a second vector orthogonal to the normal vector and
        # orthogonal to the first vector (v)
        # again, there are infinitely many solutions, varying in length
        s = solve([
            list(self.n) + [0],
            list(v) + [0],
        ])
        w = Vector(*s(1))
        return (self.p.pv(), v, w)

    def __neg__(self):
        """Return the negative HalfSpace, the normal is the negative normal"""
        return HalfSpace(-self.plane)

__all__ = ("HalfSpace",)
