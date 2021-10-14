# -*- coding: utf-8 -*-

# Example for generating a simple netgen mesh

from ngsolve import *
from netgen.geom2d import SplineGeometry

#Start geometry
geo = SplineGeometry()

#Define coordinates of interest
p2,p3,p4 = [ geo.AppendPoint(x,y) for x,y in [(5,0), (5,5), (0,5)] ]
p5,p6 = [ geo.AppendPoint(x,y) for x,y in [(5,-5), (0,-5)] ]
p7,p8,p9,p10 = [ geo.AppendPoint(x,y) for x,y in [(0,30),(0,-30),(20,-30),(20,30)] ]
#Define materials
geo.SetMaterial(1, "tumor")
geo.SetMaterial(2, "healthy")

#Define geometry:
geo.Append (["line", p6, p4],leftdomain=0, rightdomain=1, bc="lightsource")
geo.Append (["spline3", p2, p3, p4],leftdomain=1, rightdomain=2, bc="boundary")
geo.Append (["spline3", p6, p5, p2],leftdomain=1, rightdomain=2, bc="boundary")
geo.Append (["line", p7, p4], leftdomain=2, bc="boundary1")
geo.Append (["line", p6, p8], leftdomain=2, bc="boundary1")
geo.Append (["line", p8, p9], leftdomain=2, bc="boundary1")
geo.Append (["line", p9, p10], leftdomain=2, bc="boundary1")
geo.Append (["line", p10, p7], leftdomain=2, bc="boundary1")

mesh = geo.GenerateMesh(maxh=0.5)

mesh.Save("square_with_circle.vol")
