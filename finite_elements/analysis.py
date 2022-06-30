#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing objects related to finite elements analysis
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
# from matplotlib.colors import LinearSegmentedColormap
import numpy as npy
# import matplotlib.tri as mtri
import volmdlr as vm
import volmdlr.mesh as vmmesh
# import math
from scipy import sparse
# from scipy import linalg
# import time 
from dessia_common import DessiaObject
from typing import List #Tuple, TypeVar
import finite_elements.elements
from finite_elements.loads import ConstantLoad, SingleNodeLoad, MagnetLoad
from finite_elements.conditions import ContinuityCondition
from finite_elements.results import Result
from finite_elements.core import blue_red


class FiniteElementAnalysis(DessiaObject):
    """
    :param mesh: The meshing of the machine.
    :type mesh: Mesh object
    :param element_loads: The list of the loads applied to the triangluar elements.
    :type element_loads: List of ConstantLoad objects
    :param node_loads: The list of the loads applied to the nodes.
    :type node_loads: List of SingleNodeLoad objects
    :param magnet_loads: The list of the loads applied to the triangular elements.
    :type magnet_loads: List of MagnetLoad objects
    :param continuity_conditions: The list of continuity conditions applied to the nodes.
    :type continuity_conditions: List of ContinuityCondition objects
    """
    def __init__(self, mesh: vmmesh.Mesh,
                 element_loads: List[ConstantLoad],
                 node_loads: List[SingleNodeLoad],
                 magnet_loads: List[MagnetLoad],
                 continuity_conditions: List[ContinuityCondition],
                 node_boundary_conditions: List[finite_elements.conditions.NodeBoundaryCondition],
                 element_boundary_conditions: List[finite_elements.conditions.ElementBoundaryCondition]):
        self.mesh = mesh
        self.element_loads = element_loads  # current density J
        self.node_loads = node_loads 
        self.magnet_loads = magnet_loads
        self.continuity_conditions = continuity_conditions
        self.node_boundary_conditions = node_boundary_conditions
        self.element_boundary_conditions = element_boundary_conditions
        
        self.nb_loads = len(node_loads)
        
        DessiaObject.__init__(self, name='')

    def create_matrix(self):
        row_ind = []
        col_ind = []
        data = []

        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:

                data.extend(element.elementary_matrix())
                row_ind_n, col_ind_n = self.get_row_col_indices(element)

                row_ind.extend(row_ind_n)
                col_ind.extend(col_ind_n)

        positions = finite_elements.core.global_matrix_positions(dimension=element.dimension,
                                                                 nodes_number=len(self.mesh.nodes))

        for i, load in enumerate(self.node_loads):
            index = self.mesh.node_to_index[load.node]

            row_ind.extend((len(self.mesh.nodes) * self.dimension + i, positions[(index, load.dimension)]))
            col_ind.extend((positions[(index, load.dimension)], len(self.mesh.nodes) * self.dimension + i))
            data.extend((1, 1))
            
        for i, condition in enumerate(self.continuity_conditions):
            index1 = self.mesh.node_to_index[condition.node1]
            index2 = self.mesh.node_to_index[condition.node2]
            
            row_ind.extend((len(self.mesh.nodes)+len(self.node_loads)+i,
                            index1,
                            len(self.mesh.nodes)+len(self.node_loads)+i,
                            index2))
            col_ind.extend((index1,
                            len(self.mesh.nodes)+len(self.node_loads)+i,
                            index2,
                            len(self.mesh.nodes)+len(self.node_loads)+i))
            data.extend((1, 1, -condition.value, -condition.value))
            
        matrix = sparse.csr_matrix((data, (row_ind, col_ind)))
        return matrix
            
    def create_source_matrix(self):
        matrix = npy.zeros((len(self.mesh.nodes)*self.dimension+self.nb_loads+len(self.continuity_conditions), 1))
        for load in self.element_loads:
            for element in load.elements:
                indexes = [self.mesh.node_to_index[element.points[0]],
                           self.mesh.node_to_index[element.points[1]],
                           self.mesh.node_to_index[element.points[2]]]

                elementary_source_matrix = element.elementary_source_matrix(indexes)
                
                matrix[indexes[0]][0] += load.value * elementary_source_matrix[0]
                matrix[indexes[1]][0] += load.value * elementary_source_matrix[1]
                matrix[indexes[2]][0] += load.value * elementary_source_matrix[2]
                
        for i, load in enumerate(self.node_loads):
            matrix[len(self.mesh.nodes)*self.dimension + i][0] += load.value
            
        for magnet_load in self.magnet_loads:
            for linear_element in magnet_load.contour_linear_elements():
                indexes = [self.mesh.node_to_index[linear_element.points[0]],
                           self.mesh.node_to_index[linear_element.points[1]]]
                length = linear_element.length()
                dl = vm.Vector2D([-linear_element.interior_normal[1],
                                  linear_element.interior_normal[0]])
                matrix[indexes[0]][0] += magnet_load.magnetization_vector.Dot(dl) * length/2
                matrix[indexes[1]][0] += magnet_load.magnetization_vector.Dot(dl) * length/2
                
        return matrix

    @property
    def dimension(self):
        return self.mesh.elements_groups[0].elements[0].dimension

    def get_row_col_indices(self, element):

        indexes = [self.mesh.node_to_index[element.points[0]],
                   self.mesh.node_to_index[element.points[1]],
                   self.mesh.node_to_index[element.points[2]]]

        positions = finite_elements.core.global_matrix_positions(dimension=self.dimension,
                                                                 nodes_number=len(self.mesh.nodes))

        row_ind, col = [], []
        for index in indexes:
            for i in range(element.dimension):
                row_ind.extend(len(indexes)*element.dimension * [positions[(index, i+1)]])
                col.append(positions[(index, i+1)])

        col_ind = []
        for index in indexes:
            for i in range(element.dimension):
                col_ind.extend(col)

        return row_ind, col_ind

    def apply_boundary_conditions(self, rigidity_matrix, source_matrix):

        node_boundary_conditions = [(self.mesh.node_to_index[node_condition.application],
                    node_condition.dimension) for node_condition in self.node_boundary_conditions]
        node_boundary_conditions=sorted(node_boundary_conditions, key=lambda i: (i[0], i[1]), reverse=True)

        positions = finite_elements.core.global_matrix_positions(dimension=self.dimension,
                                                                 nodes_number=len(self.mesh.nodes))

        for node_condition in node_boundary_conditions:
            for i in range(rigidity_matrix.shape[1]):
                rigidity_matrix[positions[node_condition], i] = 0
            for i in range(rigidity_matrix.shape[0]):
                rigidity_matrix[i, positions[node_condition]] = 0

            for i in range(source_matrix.shape[0]):
                source_matrix[positions[node_condition]] = 0

        return rigidity_matrix, source_matrix

    def solve(self):
        """
        Solve the matix equation : F = K.X, where X is the unknown vector. \
        Each value of the vector represent the amplitude of the vector potential \
        of a node. 
        Returns a Result object.
        """
        # TEST RESULT : 
        # Between :
        # - npy.linalg.solve
        # - scipy.linalg.solve
        # - scipy.linalg.solve(assume_a='sym')
        # - scipy.sparse.linalg.spsolved
        # - scipy.pinvh than dot
        
        # scipy.sparse.linalg.spsolved is the fastest !
        # print('avant')
        K_sparse = self.create_matrix()
        F = self.create_source_matrix()

        K_sparse, F = self.apply_boundary_conditions(K_sparse, F)

        try:
            X = sparse.linalg.spsolve(K_sparse, F,
                                      permc_spec='NATURAL',
                                      use_umfpack=True)
        except sparse.linalg.MatrixRankWarning:
            print('MatricRankWarning')
            raise NotImplementedError
        X = list(X)
        # print('apres')
        return Result(self.mesh, X)
        
    def plot_elements_loads(self, ax=None):
        """ 
        Plots the mesh. The triangular elements are filled red if they are a \
        source of magnetic potential thanks to their current density.
        """
        if ax is None:
            fig, ax = plt.subplots()
        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:
                element.plot(ax=ax, color='w', fill=True)
        
        for load in self.element_loads:
            for element in load.elements:
                element.plot(ax=ax, color='r', fill=True)
        return ax
    
    def plot_elements_permeability(self, ax=None):
        """ 
        Plots the mesh with colored triangular elements depending on the \ 
        permeability of the material the element meshes. 
        """
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()
        
        color_map = ((0, 0, 1), (1, 0, 0))
        permeabilities = []
        for elements_group in self.mesh.elements_groups:
            permeabilities.append(elements_group.elements[0].mu_total)
        mu_max = max(permeabilities)
        mu_min = min(permeabilities)
        colors = []
        for elements_group in self.mesh.elements_groups:
            x = (elements_group.elements[0].mu_total - mu_min) / (mu_max - mu_min)
            color = (color_map[0][0]-(color_map[0][0]-color_map[1][0])*x, 
                      color_map[0][1]-(color_map[0][1]-color_map[1][1])*x,
                      color_map[0][2]-(color_map[0][2]-color_map[1][2])*x)
            colors.append(color)
        
        norm = mpl.colors.Normalize(vmin=mu_min, vmax=mu_max)
        sm = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)  #plt.get_cmap('jet_r')
        sm.set_array([])
        cbar = fig.colorbar(sm, ticks=npy.linspace(mu_min, mu_max, 10))
        cbar.set_label('permeability')
        for i, elements_group in enumerate(self.mesh.elements_groups):
            for element in elements_group.elements:
                element.plot(ax=ax, color=colors[i], fill=True)
                
        return ax
    
    def plot_magnet_loads(self, ax=None):
        """ 
        Plots the mesh. The triangular elements are filled red if they are a \
        source of magnetic potential thanks to their magnetization. The contour \
        of the magnetizating mesh is drawn in blue. 
        """
        if ax is None:
            fig, ax = plt.subplots()
        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:
                element.plot(ax=ax, color='w', fill=True)
            
        for magnet_load in self.magnet_loads:
            for element in magnet_load.elements:
                element.plot(ax=ax, color='r', fill=True)
                
            contour_linear_elements = magnet_load.contour_linear_elements()
            for linear_element in contour_linear_elements:
                linear_element.plot(ax=ax, color='b')
        
        print(self.magnet_loads)
        
        return ax
    
    def plot_continuity_condition(self, ax=None):
        if ax is None:
            ax = self.mesh.plot()
            
        for i, continuity_condition in enumerate(self.continuity_conditions):
            continuity_condition.node1.MPLPlot(ax=ax, color='C{}'.format(i % 10))
            continuity_condition.node2.MPLPlot(ax=ax, color='C{}'.format(i % 10))
        return ax