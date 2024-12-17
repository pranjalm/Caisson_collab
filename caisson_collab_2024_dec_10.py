# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior

dimtr = 2
mdb.Model(name='test', modelType=STANDARD_EXPLICIT)
s = mdb.models['test'].ConstrainedSketch(name='caisson', sheetSize=10.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.Line(point1=(-2.5, 5.0), point2=(-2.7, 5.0))
s.Line(point1=(-2.7, 5.0), point2=(-2.7, 0.0))
s.Line(point1=(-2.7, 0.0), point2=(-2.5, 0.0))
s.Line(point1=(-2.5, 0.0), point2=(-2.5, 5.0))
s.ConstructionLine(point1=(dimtr, -5.0), point2=(dimtr, 5.0))
s.FixedConstraint(entity=g[2])

# create part using skatch
p = mdb.models['test'].Part(name='caisson', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p = mdb.models['test'].parts['caisson']
p.BaseSolidRevolve(sketch=s, angle=180.0, flipRevolveDirection=OFF)




