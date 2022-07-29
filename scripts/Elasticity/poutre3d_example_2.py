#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 2022

@author: s.bendjebla
"""

import volmdlr.gmsh
import volmdlr.mesh as vmmesh
import finite_elements as fe
import finite_elements.elements
import finite_elements.loads
import finite_elements.analysis
import finite_elements.conditions
import time
import matplotlib.pyplot as plt

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

def find_points_with_z(nodes, z):
    indices_z = []
    for i, node in enumerate(nodes):
        if node.z == z:
            indices_z.append(i)
    return indices_z

def find_points_with_xz(nodes, x, z):
    indices = []
    for i, node in enumerate(nodes):
        if node.x == x and node.z == z:
            indices.append(i)
    return indices


# %% Mesh3D

elasticity_modulus, poisson_ratio, thickness, mass_density = 30*1e6, 0.25, 1, 2.7
# elasticity_modulus, poisson_ratio, thickness, mass_density = 70*1e6, 0.33, 1, 2700 #aluminium
# elasticity_modulus, poisson_ratio, thickness, mass_density = 210*1e9, 0.25, 1, 7860 #acier

files_path = ['../InputFiles/3D/poutre3d_0.5', '../InputFiles/3D/poutre3d_1',
              '../InputFiles/3D/poutre3d_3']

file_path = files_path[2]

gmsh = volmdlr.gmsh.Gmsh.from_file(file_path+'.msh')

mesh = gmsh.define_tetrahedron_element_mesh()

# ax= mesh.plot()
# ax.set_title(str(len(mesh.nodes))+' Nodes')
# ax.set_box_aspect([10,2,2])

# %% Finite Element Mesh2D

group_elements = []

for group in mesh.elements_groups:
    solid_elments2d = []
    for triangle in group.elements:
        solid_elments2d.append(fe.elements.ElasticityTetrahedralElement3D(
            triangle, elasticity_modulus, poisson_ratio, mass_density))

    group_elements.append(vmmesh.ElementsGroup(solid_elments2d, ''))

mesh_fe = vmmesh.Mesh(group_elements)
mesh_fe.nodes = gmsh.nodes[0]['all_nodes'] #Keep Gmsh order
mesh_fe.node_to_index = {mesh_fe.nodes[i]: i for i in range(len(mesh_fe.nodes))}

# %% Loads/Conditions

load = -1000
node_loads = []
application_indices = find_points_with_xz(nodes=mesh_fe.nodes, x=10, z=2)
for index in application_indices:
    node_loads.append(fe.loads.NodeLoad(mesh_fe.nodes[index], load, 3))


application_indices = find_points_with_x(nodes=mesh_fe.nodes, x=0)
node_boundary_conditions = []
for index in application_indices:
    node_boundary_conditions.extend([finite_elements.conditions.NodeBoundaryCondition(mesh_fe.nodes[index], 0, 1),
                                     finite_elements.conditions.NodeBoundaryCondition(mesh_fe.nodes[index], 0, 2),
                                     finite_elements.conditions.NodeBoundaryCondition(mesh_fe.nodes[index], 0, 3)])

#%% Analysis: plane_strain

# plane_strain, plane_stress = True, False

#%% Analysis: plane_stress

plane_strain, plane_stress = False, True

# %% Analysis

analysis = fe.analysis.FiniteElementAnalysis(mesh_fe, [], node_loads, [], [], node_boundary_conditions, [],
                                              plane_strain=plane_strain, plane_stress=plane_stress)

# %% Results

m = analysis.create_matrix()

results = analysis.solve()

elasticity_result = fe.results.ElasticityResults3D(results.mesh, results.result_vector, analysis.plane_strain, analysis.plane_stress)


# %% VTK files generation

elasticity_result.update_vtk_with_results(
            input_file_name = file_path+'.vtk',
            output_file_name = file_path+'_displacements'+'.vtk')

# %% Anaytic Result

b = 2
thickness = 2
I = ((b**3)*thickness)/12 #I = (b*(thickness**3) + thickness*(b**3))/12
z_displacement = (load*10**3)/(3*elasticity_modulus*I)


# %% Comparaison

application_indices = find_points_with_xz(nodes=mesh_fe.nodes, x=10, z=2)

z_displacements = []
for index in application_indices:
    z_displacements.append(elasticity_result.displacement_vectors_per_node[index].z)


print('Analytic:', z_displacement)
print('Solver:', z_displacements)
