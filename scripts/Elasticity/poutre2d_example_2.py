#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 2022

@author: s.bendjebla
"""

import volmdlr as vm
import volmdlr.mesh as mesh
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


mesh_fe = mesh.Mesh([group_phase1, group_phase2, group_phase3])
ax = mesh_fe.plot()
ax.set_title("Initial Mesh")

# %% Loads/Conditions

load = -10000000
application_index = mesh_fe.node_to_index[mesh.Node2D(9,1)]
node_loads = [fe.loads.NodeLoad(mesh_fe.nodes[application_index], load, 2)]

application_indices = find_points_with_x(nodes=mesh_fe.nodes, x=0)
node_boundary_conditions = []
for index in application_indices:
    node_boundary_conditions.extend([finite_elements.conditions.NodeBoundaryCondition(mesh_fe.nodes[index], 0, 1),
                                finite_elements.conditions.NodeBoundaryCondition(mesh_fe.nodes[index], 0, 2)])

#%% Analysis: plane_strain

plane_strain, plane_stress = True, False

#%% Analysis: plane_stress

# plane_strain, plane_stress = False, True

# %% Analysis

analysis = fe.analysis.FiniteElementAnalysis(mesh_fe, [], node_loads, [], [], node_boundary_conditions, [],
                                             plane_strain=plane_strain, plane_stress=plane_stress)

# %% Results

m = analysis.create_matrix()

results = analysis.solve()

elasticity_result = fe.results.ElasticityResults2D(results.mesh, results.result_vector, analysis.plane_strain, analysis.plane_stress)

# %% Plots

elasticity_result.plot_deformed_mesh(amplitude=10)
elasticity_result.plot_displacement_vectors_per_node(amplitude=0.2)

elasticity_result.plot_displacement_per_node_xy(amplitude=10)

# elasticity_result.plot_strain()
# elasticity_result.plot_stress()
