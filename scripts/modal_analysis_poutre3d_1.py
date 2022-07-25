#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 2022

@author: b.soumiya
"""

import volmdlr.gmsh
import volmdlr.mesh as vmmesh
import finite_elements as fe
import finite_elements.elements
import finite_elements.loads
import finite_elements.analysis
import finite_elements.conditions

# %% Mesh2D

file_path = 'poutre_3d_1.msh'

gmsh = volmdlr.gmsh.Gmsh.from_file(file_path)

mesh = gmsh.define_tetrahedron_element_mesh()

# mesh.plot()

# %% Finite Element Mesh2D

elasticity_modulus, poisson_ratio, thickness, mass_density = 210*1e9, 0.25, 1, 7.860 #acier

group_elements = []

for group in mesh.elements_groups:
    solid_elments3d = []
    for tetrahedral in group.elements:
        solid_elments3d.append(fe.elements.ElasticityTetrahedralElement3D(
            tetrahedral, elasticity_modulus, poisson_ratio, mass_density))

    group_elements.append(vmmesh.ElementsGroup(solid_elments3d, ''))

mesh = vmmesh.Mesh(group_elements)


# %% Analysis

analysis = fe.analysis.FiniteElementAnalysis(mesh, [], [], [], [], [], [])

eigvals, eigvecs = analysis.modal_analysis()
elasticity_results = []

for eigvec in eigvecs.T[0:5]:
    elasticity_results.append(fe.results.ElasticityResults2D(analysis.mesh,
                                                              eigvec))

for elasticity_result in elasticity_results[0:5]:
    # elasticity_result.plot_deformed_mesh()
    elasticity_result.plot_displacement_per_node_xy()

# %%

# import matplotlib.pyplot as plt
# b = [max(abs(eigvec)) for eigvec in eigvecs]
# plt.plot(b)

# plt.plot(eigvals)


