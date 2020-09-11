#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Milan
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as npy
import matplotlib.tri as mtri
import volmdlr as vm
import volmdlr.mesh as vmmesh
import finite_elements.electromag as els
import math
from scipy import sparse
from scipy import linalg
import time 

from dessia_common import DessiaObject
from typing import TypeVar, List, Tuple
MU = 4*math.pi*1e-7

l=0.1

p1=vm.Point2D((0,0))
p2=vm.Point2D((l,0))
p3=vm.Point2D((l,l))
p4=vm.Point2D((0,l))
p5=vm.Point2D((2*l,0))
p6=vm.Point2D((2*l,l))

tr2=vmmesh.TriangularElement([p1,p3,p4])
tr1=vmmesh.TriangularElement([p1,p2,p3])
tr3=vmmesh.TriangularElement([p2,p3,p6])
tr4=vmmesh.TriangularElement([p2,p5,p6])

load=els.ConstantLoad([tr1,tr2,tr3,tr4],50)


noad_load_1=els.SingleNodeLoad(p1,100)
noad_load_2=els.SingleNodeLoad(p2,40)
noad_load_3=els.SingleNodeLoad(p3,10)
noad_load_4=els.SingleNodeLoad(p4,10)
noad_load_5=els.SingleNodeLoad(p5,20)

Triangles=[tr1,tr2,tr3,tr4]

elements_group_1=vmmesh.ElementsGroup([tr1,tr2],'first_elements_group')
elements_group_2=vmmesh.ElementsGroup([tr3,tr4],'second_elements_group')

mesh=vmmesh.Mesh([elements_group_1,elements_group_2])




solution = els.FiniteElementAnalysis(mesh,[load],[noad_load_1,noad_load_2,noad_load_3,noad_load_4,noad_load_5],[],[])
M=solution.create_matrix()
M2=M.toarray()
print(solution.create_matrix())
result=solution.solve()
print(result.result_vector)

# #result.plot_magnetic_field_contour(),
# result.plot_potential_vector(None),
# # result.plot_magnetic_field(None,None)
# # result.plot_magnetic_field_vectors(None,0.005,None,None)

