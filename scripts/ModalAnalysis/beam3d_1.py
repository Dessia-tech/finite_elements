#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 1 2022

@author: s.bendjebla
"""

import volmdlr.gmsh_vm
import volmdlr.mesh as vmmesh
import finite_elements as fe
import finite_elements.elements
import finite_elements.loads
import finite_elements.analysis
import finite_elements.conditions
import math

# %% Mesh2D

files_path = ['../InputFiles/3D/beam3d_0.5', '../InputFiles/3D/beam3d_1',
              '../InputFiles/3D/beam3d_3']

file_path = files_path[0]

gmsh = volmdlr.gmsh_vm.GmshParser.from_file(file_path+'.msh')

mesh = gmsh.define_tetrahedron_element_mesh()

# mesh.plot()

# %% Finite Element Mesh2D

# elasticity_modulus, poisson_ratio, mass_density = 30*1e6, 0.25, 2.7
# elasticity_modulus, poisson_ratio, mass_density = 70*1e6, 0.33, 2700 #aluminium
elasticity_modulus, poisson_ratio, mass_density = 210*1e9, 0.25, 7860 #acier

group_elements = []

for group in mesh.elements_groups:
    solid_elments3d = []
    for tetrahedral in group.elements:
        solid_elments3d.append(fe.elements.ElasticityTetrahedralElement3D(
            tetrahedral, elasticity_modulus, poisson_ratio, mass_density))

    group_elements.append(vmmesh.ElementsGroup(solid_elments3d, ''))

mesh = vmmesh.Mesh(group_elements)

# %% Analysis

analysis = fe.analysis.FiniteElementAnalysis(mesh = mesh,
                                             element_loads = [],
                                             edge_loads = [],
                                             node_loads = [],
                                             magnet_loads = [],
                                             continuity_conditions = [],
                                             node_boundary_conditions = [],
                                             edge_boundary_conditions = [],
                                             element_boundary_conditions = [],
                                             plane_strain = False,
                                             plane_stress = True)

eigvals, eigvecs = analysis.modal_analysis(order='smallest', k=20)
elasticity_results = []
for eigvec in eigvecs[0:10]:
    elasticity_results.append(fe.results.ElasticityResults3D(analysis.mesh,
                                                              eigvec, analysis.plane_strain, analysis.plane_stress))
frequency = []
for eigval in eigvals:
    frequency.append((math.sqrt(abs(eigval))/(2*math.pi)))


# %% VTK files generation

# for i, elasticity_result in enumerate(elasticity_results):
#     elasticity_result.generate_vtk_file(file_path+'_mode_n°_'+str(i)+'.vtk')
