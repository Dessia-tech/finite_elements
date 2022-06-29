#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 2022

@author: s.bendjebla
"""

import volmdlr as vm
import volmdlr.mesh as mesh
import finite_elements as fe
import finite_elements.elements
import finite_elements.loads
import finite_elements.analysis

# %% Mesh2D

elasticity_modulus, poisson_ratio, thickness = 30*1e6, 0.25, 0.5

triangles = [mesh.TriangularElement2D([vm.Point2D(3,0),vm.Point2D(3,2),vm.Point2D(0,0)]),
             mesh.TriangularElement2D([vm.Point2D(0,2),vm.Point2D(0,0),vm.Point2D(3,2)])]

solid_elments2d = [fe.elements.SolidMechanicsTriangularElement2D(
    triangle, elasticity_modulus, poisson_ratio, thickness) for triangle in triangles]

# group_solid_elments2d = mesh.ElementsGroup(solid_elments2d, '')
# mesh = mesh.Mesh([group_solid_elments2d])

group_solid_elments2d = [mesh.ElementsGroup([solid_elment], '') for solid_elment in solid_elments2d]
mesh = mesh.Mesh(group_solid_elments2d)


# %%

node_loads = [] #[fe.loads.SingleNodeLoad(mesh.nodes[1], -1000)]

analysis = fe.analysis.FiniteElementAnalysis(mesh, [], node_loads, [], [])

m = analysis.create_matrix()

results = analysis.solve()

