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
import finite_elements.conditions

# %% Mesh2D

elasticity_modulus, poisson_ratio = 96, 1/3

points = [vm.Point3D(2,3,4), vm.Point3D(6,3,2), vm.Point3D(2,5,1), vm.Point3D(4,3,6)]

solid_elment3d = fe.elements.SolidMechanicsTetrahedralElement3D(
    mesh_element = mesh.TetrahedralElement(points),
    elasticity_modulus = elasticity_modulus,
    poisson_ratio = poisson_ratio)

stiffness_matrix = solid_elment3d.elementary_matrix()
