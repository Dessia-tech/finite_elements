#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 2022

@author: s.bendjebla
"""

import volmdlr.gmsh
import volmdlr.mesh as vmmesh
import finite_elements as fe
import finite_elements.elements
import finite_elements.loads
import finite_elements.analysis
import finite_elements.conditions

def find_points_with_x(nodes, x):
    indices_x = []
    for i, node in enumerate(nodes):
        if node.x == x:
            indices_x.append(i)
    return indices_x

def find_points_with_y(nodes, y):
    indices_y = []
    for i, node in enumerate(nodes):
        if node.y == y:
            indices_y.append(i)
    return indices_y

# %% Mesh2D

files_path = ['../InputFiles/2D/poutre_2d_0.8', '../InputFiles/2D/poutre_2d_0.5',
              '../InputFiles/2D/poutre_2d_0.3', '../InputFiles/2D/poutre_2d_0.18',
              '../InputFiles/2D/poutre_2d_0.1']

file_path = files_path[2]

gmsh = volmdlr.gmsh.Gmsh.from_file(file_path+'.msh')
mesh = gmsh.define_triangular_element_mesh()

# mesh.plot()

# %% Finite Element Mesh2D

elasticity_modulus, poisson_ratio, thickness, mass_density = 30*1e6, 0.25, 1, 2.7
# elasticity_modulus, poisson_ratio, thickness, mass_density = 70*1e6, 0.33, 1, 2700 #aluminium
# elasticity_modulus, poisson_ratio, thickness, mass_density = 210*1e9, 0.25, 1, 7860 #acier

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


# %% Loads/Conditions

application_indices = find_points_with_x(nodes=mesh.nodes, x=0)
# application_indices.extend(find_points_with_x(nodes=mesh.nodes, x=10))

node_boundary_conditions = []
for index in application_indices:
    node_boundary_conditions.extend([finite_elements.conditions.NodeBoundaryCondition(mesh.nodes[index], 0, 1),
                                     finite_elements.conditions.NodeBoundaryCondition(mesh.nodes[index], 0, 2)])

load = -1000

application_index = mesh.node_to_index[volmdlr.mesh.Node2D(10,1)]
node_loads = [fe.loads.NodeLoad(mesh.nodes[application_index], load, 2)]

#%% Analysis: plane_strain

# plane_strain, plane_stress = True, False

#%% Analysis: plane_stress

plane_strain, plane_stress = False, True

# %% Static Analysis

analysis = fe.analysis.FiniteElementAnalysis(mesh, [], node_loads, [], [], node_boundary_conditions, [],
                                             plane_strain=plane_strain, plane_stress=plane_stress)

# %%% Static Results

results = analysis.solve()

elasticity_result = fe.results.ElasticityResults2D(results.mesh, results.result_vector, analysis.plane_strain, analysis.plane_stress)



'''

elasticity_results_displacements[file_path[17::]] = {'nodes':len(mesh.nodes),
                                                     'displacement_vector':elasticity_result.displacement_vectors_per_node[application_index],
                                                     'time':end-start}
# y_displacements.append(elasticity_result.displacement_vectors_per_node[application_index].y)

elasticity_result.update_vtk_with_results(
            input_file_name = file_path+'.vtk',
            output_file_name = file_path+'_displacements'+'.vtk')


# %% Modal Analysis


eigvals, eigvecs = analysis.modal_analysis()
elasticity_results = []

for eigvec in eigvecs.T[0:10]:
    elasticity_results.append(fe.results.ElasticityResults2D(analysis.mesh,
                                                             eigvec,
                                                             analysis.plane_strain, analysis.plane_stress))

for elasticity_result in elasticity_results:
    # elasticity_result.plot_deformed_mesh(amplitude=50)
    elasticity_result.plot_displacement_per_node_xy(amplitude=50)

# %% VTK files generation

for i, elasticity_result in enumerate(elasticity_results):
    elasticity_result.update_vtk_with_results(
        input_file_name = file_path+'.vtk',
        output_file_name = file_path+'_mode_n°_'+str(i)+'.vtk')
'''