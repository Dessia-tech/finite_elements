#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 2022

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

files_path = ['mesh_2', 'mesh_3', 'mesh_4']

file_path = files_path[1]

gmsh = volmdlr.gmsh.Gmsh.from_file(file_path+'.msh')

mesh = gmsh.define_triangular_element_mesh()

# mesh.plot()

# %% Finite Element Mesh2D

# elasticity_modulus, poisson_ratio, thickness, mass_density = 70*1e6, 0.33, 1, 2700 #aluminium
# elasticity_modulus, poisson_ratio, thickness, mass_density = 45*1e6, 0.29, 1, 1800 #magnesium
# elasticity_modulus, poisson_ratio, thickness, mass_density = 210*1e3, 0.25, 0.5, 3.74 #acier
elasticity_modulus, poisson_ratio, thickness, mass_density = 210*1e9, 0.25, 1, 7.860 #acier
# elasticity_modulus, poisson_ratio, thickness, mass_density = 3*1e7, 0.3, 0.5, 0.3/386 #acier
# # elasticity_modulus, poisson_ratio, thickness, mass_density = 20*1e9, 0.3, 0.5, 7800 #acier
# elasticity_modulus, poisson_ratio, thickness, mass_density = 210*1e9, 0.3, 0.5, 8000 #acier

# elasticity_modulus, poisson_ratio, thickness, mass_density = 30*1e6, 0.25, 0.5, 7860
# elasticity_modulus, poisson_ratio, thickness, mass_density = 210000, 0.3, 0.5, 2.7000*1e-6 #matlab
# elasticity_modulus, poisson_ratio, thickness, mass_density = 193000, 0.3, 0.5, 8*1e-6 #Stainless steel
# elasticity_modulus, poisson_ratio, thickness, mass_density = 68900, 0.33, 0.5, 2.7*1e-6 #Aluminum 6061

group_elements = []

for group in mesh.elements_groups:
    solid_elments2d = []
    for triangle in group.elements:
        solid_elments2d.append(fe.elements.ElasticityTriangularElement2D(
            triangle, elasticity_modulus, poisson_ratio, mass_density, thickness))

    group_elements.append(vmmesh.ElementsGroup(solid_elments2d, ''))

mesh = vmmesh.Mesh(group_elements)
mesh.nodes = gmsh.nodes[0]['all_nodes'] #Keep Gmsh order
mesh.node_to_index = {mesh.nodes[i]: i for i in range(len(mesh.nodes))}


# %% Analysis

analysis = fe.analysis.FiniteElementAnalysis(mesh, [], [], [], [], [], [],
                                             plane_strain=False, plane_stress=True)

eigvals, eigvecs = analysis.modal_analysis()
elasticity_results = []

for eigvec in eigvecs.T[0:10]:
    elasticity_results.append(fe.results.ElasticityResults2D(analysis.mesh,
                                                             eigvec,
                                                             analysis.plane_strain, analysis.plane_stress))

for elasticity_result in elasticity_results:
    # elasticity_result.plot_deformed_mesh()
    elasticity_result.plot_displacement_per_node_xy()

# %% VTK files generation

for i, elasticity_result in enumerate(elasticity_results):
    elasticity_result.update_vtk_with_results(
        input_file_name = file_path+'.vtk',
        output_file_name = file_path+'_mode_n°_'+str(i)+'.vtk')
