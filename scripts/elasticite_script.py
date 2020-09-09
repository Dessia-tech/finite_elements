#!/usr/bin/env python3
# -*- coding: utf-8 -ments_gro*-
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
import finite_elements.elasticity as els
import finite_elements.core as corefe
import math
from scipy import sparse
from scipy import linalg
from finite_elements.core import steel,aluminium
import time 
from dessia_common import DessiaObject
from typing import TypeVar, List, Tuple


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
load0=els.ConstantLoad([tr1],[12,5])
load1=els.ConstantLoad([tr1,tr2],[12,5])
load2=els.ConstantLoad([tr3,tr4],[3,7])

noad_load_1=els.SingleNodeLoad(p1,[1,1])
noad_load_2=els.SingleNodeLoad(p2,[2,2])
noad_load_3=els.SingleNodeLoad(p3,[3,3])
noad_load_4=els.SingleNodeLoad(p4,[4,4])
noad_load_5=els.SingleNodeLoad(p5,[0.5,0.5])
noad_load_6=els.SingleNodeLoad(p6,[6,6])

displacement_condition_1=els.NodeDisplacement(p1,0,0)
displacement_condition_2=els.NodeDisplacement(p2,0,0) 
displacement_condition_4=els.NodeDisplacement(p4,0,0)
displacement_condition_3=els.NodeDisplacement(p3,0,0)                                       
Triangles=[tr1,tr2,tr3,tr4]

elements_group_0=vmmesh.ElementsGroup([tr1],'first_elements_group')
elements_group_1=vmmesh.ElementsGroup([tr1,tr2],'first_elements_group')
elements_group_2=vmmesh.ElementsGroup([tr3,tr4],'second_elements_group')
mesh=vmmesh.Mesh([elements_group_1,elements_group_2])
mesh.plot()

materials=corefe.Materials({elements_group_1:aluminium,elements_group_2:steel})

solution = els.FiniteElementAnalysis(mesh,materials,[load0],[noad_load_1,noad_load_2,noad_load_3,noad_load_4,noad_load_5,noad_load_6],[displacement_condition_1,displacement_condition_2],[])
M=solution.create_stiffness_matrix()
M2=M.toarray()
print(npy.linalg.matrix_rank(M2))


N=solution.create_force_matrix()
print(N)
result= solution.solve()
v=result.result_vector
print(v)



