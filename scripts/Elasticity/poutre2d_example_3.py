#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 2022

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

# %% Data

files_path = ['../InputFiles/2D/poutre_2d_0.8', '../InputFiles/2D/poutre_2d_0.5',
              '../InputFiles/2D/poutre_2d_0.3', '../InputFiles/2D/poutre_2d_0.1']

elasticity_results_displacements = {}
y_displacements = []

load = -1000

elasticity_modulus, poisson_ratio, thickness, mass_density = 30*1e6, 0.25, 1, 2.7
# elasticity_modulus, poisson_ratio, thickness, mass_density = 70*1e6, 0.33, 1, 2700 #aluminium
# elasticity_modulus, poisson_ratio, thickness, mass_density = 210*1e9, 0.25, 1, 7860 #acier

for file_path in files_path:

    # %% Mesh2D

    gmsh = volmdlr.gmsh.Gmsh.from_file(file_path+'.msh')

    mesh = gmsh.define_triangular_element_mesh()

    # mesh.plot()

    # %% Finite Element Mesh2D

    group_elements = []

    for group in mesh.elements_groups:
        solid_elments2d = []
        for triangle in group.elements:
            solid_elments2d.append(fe.elements.ElasticityTriangularElement2D(
                triangle, elasticity_modulus, poisson_ratio, mass_density, thickness))
    
        group_elements.append(vmmesh.ElementsGroup(solid_elments2d, ''))

    mesh = vmmesh.Mesh(group_elements)

    # %% Loads/Conditions

    application_index = mesh.node_to_index[volmdlr.mesh.Node2D(10,1)]
    node_loads = [fe.loads.NodeLoad(mesh.nodes[application_index], load, 2)]

    application_indices = find_points_with_x(nodes=mesh.nodes, x=0)
    node_boundary_conditions = []
    for index in application_indices:
        node_boundary_conditions.extend([finite_elements.conditions.NodeBoundaryCondition(mesh.nodes[index], 0, 1),
                                         finite_elements.conditions.NodeBoundaryCondition(mesh.nodes[index], 0, 2)])

    #%% Analysis: plane_strain

    # plane_strain, plane_stress = True, False

    #%% Analysis: plane_stress

    plane_strain, plane_stress = False, True

    # %% Analysis

    analysis = fe.analysis.FiniteElementAnalysis(mesh, [], node_loads, [], [], node_boundary_conditions, [],
                                                 plane_strain=plane_strain, plane_stress=plane_stress)

    # %% Results

    start = time.time()
    m = analysis.create_matrix()

    results = analysis.solve()
    end = time.time()

    elasticity_result = fe.results.ElasticityResults2D(results.mesh, results.result_vector, analysis.plane_strain, analysis.plane_stress)

    elasticity_results_displacements[file_path[17::]] = {'nodes':len(mesh.nodes),
                                                         'displacement_vector':elasticity_result.displacement_vectors_per_node[application_index],
                                                         'time':end-start}
    # y_displacements.append(elasticity_result.displacement_vectors_per_node[application_index].y)

    # %% Plots

    # elasticity_result.plot_deformed_mesh(amplitude=10)
    # elasticity_result.plot_displacement_vectors_per_node(amplitude=0.2)

    # elasticity_result.plot_displacement_per_node_xy(amplitude=10)

    # elasticity_result.plot_strain()
    # elasticity_result.plot_stress()


# %% Anaytic Result

b = 1
I = ((b**3)*thickness)/12 #I = (b*(thickness**3) + thickness*(b**3))/12
y_displacement = (load*10**3)/(3*elasticity_modulus*I)

# %% Comparaison

x, y1, y2 = [], [], []
for _, value in elasticity_results_displacements.items():
    x.append(value['nodes'])
    y1.append(abs(value['displacement_vector'].y-y_displacement))
    y2.append(value['time'])

# Error
fig, ax = plt.subplots()
ax.plot(x, y1, linewidth=2.0)

ax.set_xlabel('Number of Nodes')
ax.set_ylabel('Error: abs(Analytic-Solver)')

if analysis.plane_stress:
    comment = '"Plane Stress"'
else:
    comment = '"Plane Strain"'
ax.set_title(comment)

# Time
fig, ax = plt.subplots()
ax.plot(x, y2, linewidth=2.0)

ax.set_xlabel('Number of Nodes')
ax.set_ylabel('Time Computation)')

if analysis.plane_stress:
    comment = '"Plane Stress"'
else:
    comment = '"Plane Strain"'
ax.set_title(comment)

