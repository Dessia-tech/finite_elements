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


l=1

# p1=vm.Point2D((0,0))
# p2=vm.Point2D((l,0))
# p3=vm.Point2D((l,l))
# p4=vm.Point2D((0,l))
# p5=vm.Point2D((l/2,l/2))     
# p6=vm.Point2D((2*l,0))
# p6=vm.Point2D((2*l,l))

# tr2=vmmesh.TriangularElement([p1,p4,p5])
# tr1=vmmesh.TriangularElement([p1,p2,p5])
# tr3=vmmesh.TriangularElement([p2,p3,p5])
# tr4=vmmesh.TriangularElement([p3,p5,p4])
# # tr3=vmmesh.TriangularElement([p2,p3,p6])
# # tr4=vmmesh.TriangularElement([p2,p5,p6])


p1=vm.Point2D((0,0))
p2=vm.Point2D((l,0))
p3=vm.Point2D((l,l))
p4=vm.Point2D((0,l))
p5=vm.Point2D((2*l,0))
p6=vm.Point2D((2*l,l))

# tr1=vmmesh.TriangularElement([p1,p3,p5])
# tr2=vmmesh.TriangularElement([p1,p7,p5])
# tr3=vmmesh.TriangularElement([p2,p3,p6])
# tr4=vmmesh.TriangularElement([p2,p5,p6])


tr1=vmmesh.TriangularElement([p1,p2,p3])
tr2=vmmesh.TriangularElement([p1,p4,p3])
tr3=vmmesh.TriangularElement([p2,p3,p6])
tr4=vmmesh.TriangularElement([p2,p5,p6])
load0=els.ConstantLoad([tr1],[0,0])
load1=els.ConstantLoad([tr1,tr2],[12,5])
load2=els.ConstantLoad([tr3,tr4],[3,7])
boundary_load1=els.BoundaryLoad(p2,p3,vm.Vector2D([-1,0],'interior_normal'),vm.Vector2D([100,0],'load_vector'))
boundary_load2=els.BoundaryLoad(p1,p4,vm.Vector2D([-1,0],'interior_normal'),vm.Vector2D([-100,0],'load_vector'))

noad_load_1=els.SingleNodeLoad(p1,[0.1,0.1])
noad_load_2=els.SingleNodeLoad(p2,[0.1,0.2])
noad_load_3=els.SingleNodeLoad(p3,[0.3,0.3])
noad_load_4=els.SingleNodeLoad(p4,[0.,0.4])
noad_load_5=els.SingleNodeLoad(p5,[0.5,0.5])
noad_load_6=els.SingleNodeLoad(p6,[0.6,0.6])
[noad_load_1,noad_load_2,noad_load_3,noad_load_4,noad_load_5,noad_load_6]

displacement_condition_1=els.NodeDisplacement(p1,0,0)
displacement_condition_2=els.NodeDisplacement(p2,0,0) 
displacement_condition_4=els.NodeDisplacement(p4,0,0)
displacement_condition_3=els.NodeDisplacement(p3,0,0)
displacement_condition_5=els.NodeDisplacement(p5,0,0)
displacement_condition_6=els.NodeDisplacement(p6,0,0)                                           


elements_group_0=vmmesh.ElementsGroup([tr1],'first_elements_group')
elements_group_1=vmmesh.ElementsGroup([tr1,tr2],'first_elements_group')
elements_group_2=vmmesh.ElementsGroup([tr3,tr4],'second_elements_group')
mesh=vmmesh.Mesh([elements_group_1,elements_group_2])
ax=mesh.plot()

materials=corefe.Materials({elements_group_1:steel,elements_group_2:steel})

solution = els.FiniteElementAnalysis(mesh,materials,[],[boundary_load1,boundary_load2],[],[displacement_condition_1,displacement_condition_4])
M=solution.create_stiffness_matrix()
N=solution.create_force_matrix()
print(N)
M2=M.toarray()
# print(sparse.linalg.eigs(M))
# L=[420067.14909799      , 421085.77640813    ,
#        151366.43137664       ,  82840.23669258,
#         82840.23669258,  59298.56395972]
# for k in range(len(L)):
#     print(math.sqrt(L[k])/(2*math.pi))



result= solution.solve()
v=result.result_vector
print(v)
node_deformation=result.node_to_displacement()
print(node_deformation)
mesh.plot_displaced_mesh(node_deformation,ax)



