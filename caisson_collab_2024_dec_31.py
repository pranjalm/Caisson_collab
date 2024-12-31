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

# create a part using extrusion HOMEWORK (10 meter cube soil part)

length = 10.0
# 创建草图，草图将会在X-Y平面
s = mdb.models['test'].ConstrainedSketch(name='soil', sheetSize=20.0)
# 在草图中绘制一个矩形，矩形的边长为10m
s.rectangle(point1=(0.0, 0.0), point2=(length, length))
# 创建立方体的零件
p = mdb.models['test'].Part(name='soil', dimensionality=THREE_D, type=DEFORMABLE_BODY)
# 使用草图来创建立方体的初步形状
p = mdb.models['test'].parts['soil']
p.BaseSolidExtrude(sketch=s, depth=length)

# created the steel material as elastic

dict_mat = {'soil':[2000, 30000000, 0.45], 'steel':[7850, 210000000000, 0.3]} 
for k,v in dict_mat.items():
    mdb.models['test'].Material(name=k) # define material name
    mdb.models['test'].materials[k].Density(table=((v[0], ), ))  # assign density
    mdb.models['test'].materials[k].Elastic(table=((v[1], v[2]), )) # assign E and mu
    mdb.models['test'].HomogeneousSolidSection(name=k+'_section', material=k, thickness=None) # create section for parts

p = mdb.models['test'].parts['caisson']
p.Set(cells= p.cells.getSequenceFromMask(('[#1 ]', ), ), name='main')
p.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
    region= p.sets['main'], sectionName='steel_section', thicknessAssignment=FROM_SECTION)

# create a section and set for 10 meter cube soil part
p = mdb.models['test'].parts['soil']
p.Set(cells= p.cells.getSequenceFromMask(('[#1 ]', ), ), name='main')
p.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
    region= p.sets['main'], sectionName='soil_section', thicknessAssignment=FROM_SECTION)

tp = 2.0
mdb.models['test'].ExplicitDynamicsStep(name='exp', timePeriod=tp, previous='Initial', improvedDtMethod=ON)
mdb.models['test'].fieldOutputRequests['F-Output-1'].setValues(variables=('A','U', 'V', 'S'), timeInterval=tp/100, timeMarks=ON)
mdb.models['test'].historyOutputRequests['H-Output-1'].setValues(variables=('ALLAE', ), timeInterval=tp/100)

# assembaly of the parts

a = mdb.models['test'].rootAssembly
for i in mdb.models['test'].parts.keys():
    p = mdb.models['test'].parts[i]
    a.Instance(name=i+'-1', part=p, dependent=ON)
a.translate(instanceList=('caisson-1', ), vector=(3, 10, 0))

# Use a loop to create mesh to the parts as we did in the assembly
#########################################
#########################################
#########################################

p = mdb.models['test'].parts['caisson']
p.seedPart(size=0.5, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['test'].parts['caisson']
p.generateMesh()

p = mdb.models['test'].parts['soil']
p.seedPart(size=0.5, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['test'].parts['soil']
p.generateMesh()

