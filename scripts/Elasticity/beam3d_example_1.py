#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 05 2022

@author: s.bendjebla
"""

import volmdlr as vm
import volmdlr.mesh as mesh
import finite_elements as fe
import finite_elements.elements
import finite_elements.loads
import finite_elements.analysis
import finite_elements.conditions

# %% Mesh3D

# elasticity_modulus, poisson_ratio, thickness, mass_density = 70*1e6, 0.33, 1, 2700 #aluminium
elasticity_modulus, poisson_ratio, mass_density = 210*1e9, 0.25, 7860 #acier

points = [vm.Point3D(2,3,4), vm.Point3D(6,3,2), vm.Point3D(2,5,1), vm.Point3D(4,3,6)]

solid_elments3d = [fe.elements.ElasticityTetrahedralElement3D(
    mesh_element = mesh.TetrahedralElement(points),
    elasticity_modulus = elasticity_modulus,
    poisson_ratio = poisson_ratio,
    mass_density = mass_density)]

group_solid_elments3d = mesh.ElementsGroup(solid_elments3d, '')
mesh = mesh.Mesh([group_solid_elments3d])

stiffness_matrix = solid_elments3d[0].elementary_matrix(plane_strain=True, plane_stress=False)


# %% Analysis

analysis = fe.analysis.FiniteElementAnalysis(mesh = mesh,
                                             element_loads = [],
                                             edge_loads = [],
                                             node_loads = [] ,
                                             magnet_loads = [],
                                             continuity_conditions = [],
                                             node_boundary_conditions = [],
                                             edge_boundary_conditions = [],
                                             element_boundary_conditions = [],
                                             plane_strain = True,
                                             plane_stress = False)

m = analysis.create_matrix()
