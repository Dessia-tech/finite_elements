#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 2022

@author: b.soumiya
"""

import volmdlr as vm
import volmdlr.mesh as mesh
import finite_elements as fe
import finite_elements.elements
import finite_elements.loads
import finite_elements.analysis
import finite_elements.conditions

# %% Mesh2D

phase1, phase2, phase3 = 3, 3, 3
# elasticity_modulus, poisson_ratio, thickness, mass_density = 70*1e6, 0.33, 1, 2700 #aluminium
elasticity_modulus, poisson_ratio, thickness, mass_density = 210*1e9, 0.25, 1, 7860 #acier


elements_phase1 = []
for i in range(3):
    element1 = mesh.TriangularElement2D([mesh.Node2D(i,0), mesh.Node2D(i+1,0), mesh.Node2D(i,1)])
    element2 = mesh.TriangularElement2D([mesh.Node2D(i+1,1), mesh.Node2D(i+1,0), mesh.Node2D(i,1)])
    elements_phase1.extend([element1, element2])

solid_elments2d = [fe.elements.ElasticityTriangularElement2D(
    element, elasticity_modulus, poisson_ratio, mass_density, thickness) for element in elements_phase1]
group_phase1 = mesh.ElementsGroup(solid_elments2d, 'phase1')


elements_phase2 = []
for i in range(phase1, phase1+phase2):
    element1 = mesh.TriangularElement2D([mesh.Node2D(i,0), mesh.Node2D(i+1,0), mesh.Node2D(i,1)])
    element2 = mesh.TriangularElement2D([mesh.Node2D(i+1,1), mesh.Node2D(i+1,0), mesh.Node2D(i,1)])
    elements_phase2.extend([element1, element2])

solid_elments2d = [fe.elements.ElasticityTriangularElement2D(
    element, elasticity_modulus, poisson_ratio, mass_density, thickness) for element in elements_phase2]
group_phase2 = mesh.ElementsGroup(solid_elments2d, 'phase2')


elements_phase3 = []
for i in range(phase1+phase2, phase1+phase2+phase3):
    element1 = mesh.TriangularElement2D([mesh.Node2D(i,0), mesh.Node2D(i+1,0), mesh.Node2D(i,1)])
    element2 = mesh.TriangularElement2D([mesh.Node2D(i+1,1), mesh.Node2D(i+1,0), mesh.Node2D(i,1)])
    elements_phase3.extend([element1, element2])

solid_elments2d = [fe.elements.ElasticityTriangularElement2D(
    element, elasticity_modulus, poisson_ratio, mass_density, thickness) for element in elements_phase3]
group_phase3 = mesh.ElementsGroup(solid_elments2d, 'phase3')


mesh = mesh.Mesh([group_phase1, group_phase2, group_phase3])

# %% Analysis

analysis = fe.analysis.FiniteElementAnalysis(mesh, [], [], [], [], [], [],
                                             plane_strain=False, plane_stress=True)

eigvals, eigvecs = analysis.modal_analysis()
elasticity_results = []

for eigvec in eigvecs.T:
    elasticity_results.append(fe.results.ElasticityResults2D(analysis.mesh,
                                                             eigvec,
                                                             analysis.plane_strain, analysis.plane_stress))

for elasticity_result in elasticity_results[0:10]:
    # elasticity_result.plot_deformed_mesh()
    elasticity_result.plot_displacement_per_node_xy()

# %%

# import matplotlib.pyplot as plt
# b = [max(abs(eigvec)) for eigvec in eigvecs]
# plt.plot(b)

# plt.plot(eigvals)

