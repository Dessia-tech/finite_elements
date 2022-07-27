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

file_path = 'poutre_0.5.msh'

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

analysis = fe.analysis.FiniteElementAnalysis(mesh, [], [], [], [], [], [],
                                             plane_strain=False, plane_stress=True)


eigvals, eigvecs = analysis.modal_analysis()
elasticity_results = []

for eigvec in eigvecs.T[0:10]:
    elasticity_results.append(fe.results.ElasticityResults3D(analysis.mesh,
                                                              eigvec, analysis.plane_strain, analysis.plane_stress))

import matplotlib.pyplot as plt

for elasticity_result in elasticity_results[0:5]:
    # elasticity_result.plot_deformed_mesh()
    # elasticity_result.plot_displacement_per_node_xy()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1,projection='3d')
    ax.set_box_aspect([10,2,2])
    elasticity_result.deformed_mesh.plot(ax)


# %%

# import matplotlib.pyplot as plt
# b = [max(abs(eigvec)) for eigvec in eigvecs]
# plt.plot(b)

# plt.plot(eigvals)

# %% get initial nodes

nodes = gmsh.nodes[0]['all_nodes']
i=0
mesh_nodes = mesh.nodes
i=i+1

lines = []
for node in nodes:
    index = mesh.node_to_index[node]
    lines.append(str(elasticity_results[i].displacement_vectors_per_node[index].norm()))

with open('poutre_displacement' + '.vtk', 'w', encoding="utf-8") as f:
    for line in lines:
        f.write(line)
        f.write('\n')
f.close()


lines = []
for node in nodes:
    index = mesh.node_to_index[node]
    n=elasticity_results[i].deformed_nodes[index]
    lines.append(str(n.x)+' '+str(n.y)+' '+str(n.z))

with open('poutre_nodes' + '.vtk', 'w', encoding="utf-8") as f:
    for line in lines:
        f.write(line)
        f.write('\n')
f.close()


# %%

lines = []
for displacement in elasticity_results[i].displacement_vectors_per_node:
    lines.append(str(displacement.norm()))

with open('poutre_displacement' + '.vtk', 'w', encoding="utf-8") as f:
    for line in lines:
        f.write(line)
        f.write('\n')
f.close()

# %%


