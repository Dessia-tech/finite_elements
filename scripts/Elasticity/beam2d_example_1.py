#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 2022

@author: s.bendjebla
"""

import volmdlr.mesh as mesh
import finite_elements as fe
import finite_elements.elements
import finite_elements.loads
import finite_elements.analysis
import finite_elements.conditions

#%% Mesh2D

elasticity_modulus, poisson_ratio, thickness, mass_density = 30*1e6, 0.25, 0.5, 2.7

triangles = [mesh.TriangularElement2D([mesh.Node2D(3,0),mesh.Node2D(3,2),mesh.Node2D(0,0)]),
             mesh.TriangularElement2D([mesh.Node2D(0,2),mesh.Node2D(0,0),mesh.Node2D(3,2)])]

solid_elments2d = [fe.elements.ElasticityTriangularElement2D(
    triangle, elasticity_modulus, poisson_ratio, mass_density, thickness) for triangle in triangles]

group_solid_elments2d = mesh.ElementsGroup(solid_elments2d, '')
mesh_fe = mesh.Mesh([group_solid_elments2d])
ax = mesh_fe.plot()
ax.set_title("Initial Mesh")

# %% Loads/Conditions

load = -1000 #-1000000
node_loads = [fe.loads.NodeLoad(mesh_fe.nodes[mesh_fe.node_to_index[mesh.Node2D(3,2)]], load, 2)] #1000

node_boundary_conditions = [finite_elements.conditions.NodeBoundaryCondition(mesh.Node2D(3,0), 0, 2),
                            finite_elements.conditions.NodeBoundaryCondition(mesh.Node2D(0,2), 0, 1),
                            finite_elements.conditions.NodeBoundaryCondition(mesh.Node2D(0,2), 0, 2),
                            finite_elements.conditions.NodeBoundaryCondition(mesh.Node2D(0,0), 0, 1),
                            finite_elements.conditions.NodeBoundaryCondition(mesh.Node2D(0,0), 0, 2)]


#%% Analysis: plane_strain

# plane_strain, plane_stress = True, False

#%% Analysis: plane_stress

plane_strain, plane_stress = False, True

# %% Analysis

analysis = fe.analysis.FiniteElementAnalysis(mesh = mesh_fe,
                                             element_loads = [],
                                             edge_loads = [],
                                             node_loads = node_loads ,
                                             magnet_loads = [],
                                             continuity_conditions = [],
                                             node_boundary_conditions = node_boundary_conditions,
                                             edge_boundary_conditions = [],
                                             element_boundary_conditions = [],
                                             plane_strain = plane_strain,
                                             plane_stress = plane_stress)


# %% Results

m = analysis.create_matrix()

results = analysis.solve()

elasticity_result = fe.results.ElasticityResults2D(results.mesh, results.result_vector, analysis.plane_strain, analysis.plane_stress)

# %% Plots

elasticity_result.plot_deformed_mesh(amplitude=2e3)
elasticity_result.plot_displacement_vectors_per_node(amplitude=0.1)

elasticity_result.plot_displacement_per_node_xy(amplitude=2e3)
elasticity_result.plot_energy(amplitude=2e3)

# elasticity_result.plot_strain()
# elasticity_result.plot_stress()
