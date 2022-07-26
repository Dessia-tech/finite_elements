#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 2022

@author: b.soumiya
"""

import volmdlr.gmsh
import volmdlr as vm
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

file_path = 'mesh_4.msh'

gmsh = volmdlr.gmsh.Gmsh.from_file(file_path)

mesh = gmsh.define_triangular_element_mesh()

# mesh.plot()

# %% Finite Element Mesh2D

elasticity_modulus, poisson_ratio, thickness, mass_density = 30*1e6, 0.25, 0.5, 2.7
# elasticity_modulus, poisson_ratio, thickness, mass_density = 70*1e3, 0.25, 0.5, 2.7

group_elements = []

for group in mesh.elements_groups:
    solid_elments2d = []
    for triangle in group.elements:
        solid_elments2d.append(fe.elements.ElasticityTriangularElement2D(
            triangle, elasticity_modulus, poisson_ratio, mass_density, thickness))

    group_elements.append(vmmesh.ElementsGroup(solid_elments2d, ''))

mesh = vmmesh.Mesh(group_elements)

# %% Loads/Conditions

load = -1000
application_index = mesh.node_to_index[vm.Point2D(10,1)]
node_loads = [fe.loads.NodeLoad(mesh.nodes[application_index], load, 2)]

application_indices = find_points_with_x(nodes=mesh.nodes, x=0)
node_boundary_conditions = []
for index in application_indices:
    node_boundary_conditions.extend([finite_elements.conditions.NodeBoundaryCondition(mesh.nodes[index], 0, 1),
                                finite_elements.conditions.NodeBoundaryCondition(mesh.nodes[index], 0, 2)])

#%% Analysis: plane_strain

analysis = fe.analysis.FiniteElementAnalysis(mesh, [], node_loads, [], [], node_boundary_conditions, [],
                                             plane_strain=True, plane_stress=False)

m = analysis.create_matrix()

results = analysis.solve()

elasticity_result = fe.results.ElasticityResults2D(results.mesh, results.result_vector, analysis.plane_strain, analysis.plane_stress)

elasticity_result.plot_deformed_mesh()
elasticity_result.plot_displacement_vectors_per_node()

# elasticity_result.plot_strain()
# elasticity_result.plot_stress()

# elasticity_result.plot_displacement_per_node_x()
# elasticity_result.plot_displacement_per_node_y()
elasticity_result.plot_displacement_per_node_xy()

#%% Analysis: plane_stress

analysis = fe.analysis.FiniteElementAnalysis(mesh, [], node_loads, [], [], node_boundary_conditions, [],
                                             plane_strain=False, plane_stress=True)

m = analysis.create_matrix()

results = analysis.solve()

elasticity_result = fe.results.ElasticityResults2D(results.mesh, results.result_vector, analysis.plane_strain, analysis.plane_stress)

elasticity_result.plot_deformed_mesh()
elasticity_result.plot_displacement_vectors_per_node()

# elasticity_result.plot_strain()
# elasticity_result.plot_stress()

# elasticity_result.plot_displacement_per_node_x()
# elasticity_result.plot_displacement_per_node_y()
elasticity_result.plot_displacement_per_node_xy()


# %%

b = 1
I = ((b**3)*thickness)/12 # I = (b*(thickness**3) + thickness*(b**3))/12
fleche_y = (load*10**3)/(3*elasticity_modulus*I)

elasticity_result.displacement_vectors_per_node[283]