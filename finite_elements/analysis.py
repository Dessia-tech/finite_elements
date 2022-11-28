#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing objects related to finite elements analysis
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as npy
import volmdlr as vm
import volmdlr.mesh as vmmesh
from scipy import sparse
from dessia_common import DessiaObject
from typing import List
import finite_elements.elements
from finite_elements.loads import ElementsLoad, EdgeLoad, NodeLoad, MagnetLoad
from finite_elements.conditions import ContinuityCondition
from finite_elements.results import Result
from finite_elements.core import blue_red


def node_boundary_conditions_to_dict(node_boundary_conditions):
    node_boundary_conditions_dict = {}
    for node_bc in node_boundary_conditions:
        try:
            node_boundary_conditions_dict[(node_bc.application,
                                           node_bc.dimension)] = + node_bc.value
        except KeyError:
            node_boundary_conditions_dict[(node_bc.application,
                                           node_bc.dimension)] = node_bc.value
    return node_boundary_conditions_dict


def node_boundary_from_dict(node_boundary_conditions_dict):
    node_boundary_conditions = []
    for key, value in node_boundary_conditions_dict.items():
        node_boundary_conditions.append(
            finite_elements.conditions.NodeBoundaryCondition(
                application=key[0],
                value=value,
                dimension=key[1]))
    return node_boundary_conditions


def node_loads_to_dict(node_loads):
    node_loads_dict = {}
    for node_load in node_loads:
        try:
            node_loads_dict[(node_load.node,
                             node_load.dimension)] = + node_load.value
        except KeyError:
            node_loads_dict[(node_load.node,
                             node_load.dimension)] = node_load.value
    return node_loads_dict


def node_loads_from_dict(node_loads_dict):
    node_loads = []
    for key, value in node_loads_dict.items():
        node_loads.append(NodeLoad(
            node=key[0],
            value=value,
            dimension=key[1]))
    return node_loads


class FiniteElements(DessiaObject):
    def __init__(self, mesh: vmmesh.Mesh,
                 element_loads: List[ElementsLoad],
                 edge_loads: List[EdgeLoad],
                 node_loads: List[NodeLoad],
                 magnet_loads: List[MagnetLoad],
                 continuity_conditions: List[ContinuityCondition],
                 node_boundary_conditions: List[finite_elements.conditions.NodeBoundaryCondition],
                 edge_boundary_conditions: List[finite_elements.conditions.EdgeBoundaryCondition],
                 element_boundary_conditions: List[finite_elements.conditions.ElementBoundaryCondition],
                 plane_strain: bool = None,
                 plane_stress: bool = None):
        self.mesh = mesh
        self.element_loads = element_loads  # current density J
        self.edge_loads = edge_loads
        self.node_loads = node_loads
        self.magnet_loads = magnet_loads
        self.continuity_conditions = continuity_conditions
        self.node_boundary_conditions = node_boundary_conditions
        self.edge_boundary_conditions = edge_boundary_conditions
        self.element_boundary_conditions = element_boundary_conditions
        self.plane_strain = plane_strain
        self.plane_stress = plane_stress

        self._boundary_conditions = None
        self._node_loads = None
        self._positions = None

        DessiaObject.__init__(self, name='')

    def c_matrix_continuity_conditions(self):
        row_ind, col_ind, data = [], [], []
        for i, condition in enumerate(self.continuity_conditions):
            data.extend(condition.c_matrix())

            index1 = self.mesh.node_to_index[condition.node1]
            index2 = self.mesh.node_to_index[condition.node2]

            row_ind.extend((len(self.mesh.nodes) + len(self.node_loads) + i,
                            index1,
                            len(self.mesh.nodes) + len(self.node_loads) + i,
                            index2))
            col_ind.extend((index1,
                            len(self.mesh.nodes) + len(self.node_loads) + i,
                            index2,
                            len(self.mesh.nodes) + len(self.node_loads) + i))

        return data, row_ind, col_ind

    def boundary_conditions_element_to_node(self):
        element_to_node_boundary_conditions = []
        for i, element_condition in enumerate(self.element_boundary_conditions):

            matrix_factors = element_condition.application.element_to_node_factors()

            for p, point in enumerate(element_condition.application.points):
                element_to_node_boundary_conditions.append(
                    finite_elements.conditions.NodeBoundaryCondition(
                        application=point,
                        value=element_condition.value * matrix_factors[p],
                        dimension=element_condition.dimension))
        return element_to_node_boundary_conditions

    def boundary_conditions_edge_to_node(self):
        edge_to_node_boundary_conditions = []
        for i, edge_condition in enumerate(self.edge_boundary_conditions):
            for point in [edge_condition.application.start,
                          edge_condition.application.end]:
                edge_to_node_boundary_conditions.append(
                    finite_elements.conditions.NodeBoundaryCondition(
                        application=point,
                        value=edge_condition.value * 0.5,
                        dimension=edge_condition.dimension))
        return edge_to_node_boundary_conditions

    def c_matrix_boundary_conditions(self):
        if self._boundary_conditions:
            node_boundary_conditions = self._boundary_conditions
        else:
            node_boundary_conditions = self.node_boundary_conditions[:]
            # element to node
            node_boundary_conditions.extend(self.boundary_conditions_element_to_node())
            # edge to node
            node_boundary_conditions.extend(self.boundary_conditions_edge_to_node())

            # node_bc to_dict
            node_boundary_conditions_dict = node_boundary_conditions_to_dict(node_boundary_conditions)

            # node_bc dict from_dict
            node_boundary_conditions = node_boundary_from_dict(node_boundary_conditions_dict)
            self._boundary_conditions = node_boundary_conditions

        # c_matrix data
        # positions = finite_elements.core.global_matrix_positions(dimension=self.dimension,
        #                                                          nodes_number=len(self.mesh.nodes))
        positions = self.positions
        row_ind, col_ind, data = [], [], []
        for i, node_condition in enumerate(node_boundary_conditions):
            data.extend(node_condition.c_matrix())
            pos = positions[(self.mesh.node_to_index[node_condition.application],
                             node_condition.dimension)]
            row_ind.extend((len(self.mesh.nodes) * self.dimension + i, pos))
            col_ind.extend((pos, len(self.mesh.nodes) * self.dimension + i))

        return data, row_ind, col_ind

    @property
    def elements_name(self):
        return self.mesh.elements_groups[0].elements[0].__class__.__name__

    def elements_permeability(self):

        permeabilities = [elements_group.elements[0].mu_total for elements_group in self.mesh.elements_groups]

        return permeabilities

    def k_matrix(self):
        row_ind, col_ind, data = [], [], []
        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:
                if isinstance(element, finite_elements.elements.MagneticElement2D):
                    data.extend(element.elementary_matrix())
                elif isinstance(element, finite_elements.elements.ElasticityElement):
                    data.extend(
                        element.elementary_matrix(
                            plane_strain=self.plane_strain,
                            plane_stress=self.plane_stress))

                row_ind_n, col_ind_n = self.get_row_col_indices(element)
                row_ind.extend(row_ind_n)
                col_ind.extend(col_ind_n)
        return data, row_ind, col_ind

    def m_matrix(self):
        row_ind, col_ind, data = [], [], []
        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:
                data.extend(element.elementary_mass_matrix())
                row_ind_n, col_ind_n = self.get_row_col_indices(element)
                row_ind.extend(row_ind_n)
                col_ind.extend(col_ind_n)
        return data, row_ind, col_ind

    def loads_element_to_node(self):
        element_to_node_loads = []
        for i, elements_load in enumerate(self.element_loads):
            loads_per_element = elements_load.value_per_element

            for j, element in enumerate(elements_load.elements):

                matrix_factors = element.element_to_node_factors()

                for p, point in enumerate(element.points):
                    element_to_node_loads.append(NodeLoad(
                        node=point,
                        value=loads_per_element[j] * matrix_factors[p],
                        dimension=elements_load.dimension))

        return element_to_node_loads

    def loads_edge_to_node(self):
        edge_to_node_loads = []
        for i, edge_load in enumerate(self.edge_loads):
            for point in [edge_load.edge.start,
                          edge_load.edge.end]:
                edge_to_node_loads.append(NodeLoad(
                        node=point,
                        value=edge_load.value * 0.5,
                        dimension=edge_load.dimension))
        return edge_to_node_loads

    @property
    def positions(self):
        if not self._positions:
            self._positions = finite_elements.core.global_matrix_positions(dimension=self.dimension,
                                                                           nodes_number=len(self.mesh.nodes))

        return self._positions

    def source_c_matrix_loads(self):
        node_loads = self.node_loads[:]
        # element to node
        node_loads.extend(self.loads_element_to_node())
        # edge to node
        node_loads.extend(self.loads_edge_to_node())

        # node_bc to_dict
        node_loads_dict = node_loads_to_dict(node_loads)

        # node_bc dict from_dict
        node_loads = node_loads_from_dict(node_loads_dict)

        if not self._node_loads:
            self._node_loads = node_loads

        # source_c_matrix data
        data, row_ind = [], []
        # positions = finite_elements.core.global_matrix_positions(dimension=self.dimension,
        #                                                          nodes_number=len(self.mesh.nodes))
        positions = self.positions
        for i, load in enumerate(node_loads):
            data.append(load.source_c_matrix())
            row_ind.append(positions[(self.mesh.node_to_index[load.node], load.dimension)])

        return data, row_ind

    # def source_c_matrix_elements_loads(self):
    #     data, row_ind = [], []
    #     for load in self.element_loads:
    #         for element in load.elements:
    #             indexes = [self.mesh.node_to_index[point] for point in element.points]

    #             elementary_source_matrix = element.element_to_node_factors()

    #             data.append(load.value * elementary_source_matrix[0])
    #             data.append(load.value * elementary_source_matrix[1])
    #             data.append(load.value * elementary_source_matrix[2])

    #             row_ind.append(indexes[0])
    #             row_ind.append(indexes[1])
    #             row_ind.append(indexes[2])

    #     return data, row_ind

    def source_c_matrix_boundary_conditions(self):
        if self._boundary_conditions:
            node_boundary_conditions = self._boundary_conditions
        else:
            node_boundary_conditions = self.node_boundary_conditions[:]
            # element to node
            node_boundary_conditions.extend(self.boundary_conditions_element_to_node())
            # edge to node
            node_boundary_conditions.extend(self.boundary_conditions_edge_to_node())

            # node_bc to_dict
            node_boundary_conditions_dict = node_boundary_conditions_to_dict(node_boundary_conditions)

            # node_bc dict from_dict
            node_boundary_conditions = node_boundary_from_dict(node_boundary_conditions_dict)

            self._boundary_conditions = node_boundary_conditions

        # source_c_matrix data
        data, row_ind = [], []

        for i, node_condition in enumerate(node_boundary_conditions):
            data.append(node_condition.source_c_matrix())
            row_ind.append(len(self.mesh.nodes) * self.dimension + i)

        return data, row_ind

    # def source_c_matrix_node_loads(self):
    #     data, row_ind = [], []
    # #   positions = finite_elements.core.global_matrix_positions(dimension=self.dimension,
    # #                                                            nodes_number=len(self.mesh.nodes))
    #     positions = self.positions
    #     for i, load in enumerate(self.node_loads):
    #         data.append(load.source_c_matrix())
    #         row_ind.append(positions[(self.mesh.node_to_index[load.node], load.dimension)])

    #     return data, row_ind

    def source_c_matrix_magnet_loads(self):
        data, row_ind = [], []
        for magnet_load in self.magnet_loads:
            for linear_element in magnet_load.contour_linear_elements():
                indexes = [self.mesh.node_to_index[linear_element.points[0]],
                           self.mesh.node_to_index[linear_element.points[1]]]
                length = linear_element.length()
                dl = vm.Vector2D([-linear_element.interior_normal[1],
                                  linear_element.interior_normal[0]])
                data.append(magnet_load.magnetization_vector.Dot(dl) * length / 2)
                data.append(magnet_load.magnetization_vector.Dot(dl) * length / 2)
                row_ind.append(indexes[0])
                row_ind.append(indexes[1])

        return data, row_ind


class FiniteElementAnalysis(FiniteElements):
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

    # def __init__(self, mesh: vmmesh.Mesh,
    #              element_loads: List[ElementsLoad],
    #              node_loads: List[NodeLoad],
    #              magnet_loads: List[MagnetLoad],
    #              continuity_conditions: List[ContinuityCondition],
    #              node_boundary_conditions: List[finite_elements.conditions.NodeBoundaryCondition],
    #              element_boundary_conditions: List[finite_elements.conditions.ElementBoundaryCondition]):
    #     self.mesh = mesh
    #     self.element_loads = element_loads  # current density J
    #     self.node_loads = node_loads
    #     self.magnet_loads = magnet_loads
    #     self.continuity_conditions = continuity_conditions
    #     self.node_boundary_conditions = node_boundary_conditions
    #     self.element_boundary_conditions = element_boundary_conditions

    #     self.nb_loads = len(node_loads)

    #     FiniteElements.__init__(self, mesh, element_loads, node_loads,
    #                             magnet_loads, continuity_conditions,
    #                             node_boundary_conditions, element_boundary_conditions)

    def create_matrix(self):
        row_ind = []
        col_ind = []
        data = []

        # # global K
        # k_matrix = self.k_matrix()
        # data.extend(k_matrix[0])
        # row_ind.extend(k_matrix[1])
        # col_ind.extend(k_matrix[2])

        # # continuity_conditions
        # c_matrix_continuity_conditions = self.c_matrix_continuity_conditions()
        # data.extend(c_matrix_continuity_conditions[0])
        # row_ind.extend(c_matrix_continuity_conditions[1])
        # col_ind.extend(c_matrix_continuity_conditions[2])

        # # boundary_conditions
        # c_matrix_boundary_conditions = self.c_matrix_boundary_conditions()
        # data.extend(c_matrix_boundary_conditions[0])
        # row_ind.extend(c_matrix_boundary_conditions[1])
        # col_ind.extend(c_matrix_boundary_conditions[2])

        method_names = ['k_matrix', 'c_matrix_continuity_conditions',
                        'c_matrix_boundary_conditions']

        for method_name in method_names:
            if hasattr(self, method_name):
                result = getattr(self, method_name)()
                data.extend(result[0])
                row_ind.extend(result[1])
                col_ind.extend(result[2])

            else:
                raise NotImplementedError(
                    f'Class {self.__class__.__name__} does not implement {method_name}')

        matrix = sparse.csr_matrix((data, (row_ind, col_ind)))

        return matrix

    def create_source_matrix(self):
        matrix = npy.zeros((self.get_source_matrix_length(), 1))

        # # elements_loads
        # source_c_matrix_elements_loads = self.source_c_matrix_elements_loads()
        # for i, d in enumerate(source_c_matrix_elements_loads[0]):
        #     matrix[source_c_matrix_elements_loads[1][i]][0] += d

        # # node_loads
        # source_c_matrix_node_loads = self.source_c_matrix_node_loads()
        # for i, d in enumerate(source_c_matrix_node_loads[0]):
        #     matrix[source_c_matrix_node_loads[1][i]][0] += d

        # # magnet_loads
        # source_c_matrix_magnet_loads = self.source_c_matrix_magnet_loads()
        # for i, d in enumerate(source_c_matrix_magnet_loads[0]):
        #     matrix[source_c_matrix_magnet_loads[1][i]][0] += d

        # # boundary_conditions
        # source_c_matrix_node_boundary_conditions = self.source_c_matrix_node_boundary_conditions()
        # for i, d in enumerate(source_c_matrix_node_boundary_conditions[0]):
        #     matrix[source_c_matrix_node_boundary_conditions[1][i]][0] += d

        method_names = ['source_c_matrix_loads',
                        'source_c_matrix_magnet_loads', 'source_c_matrix_boundary_conditions']

        for method_name in method_names:
            if hasattr(self, method_name):
                result = getattr(self, method_name)()
                for i, d in enumerate(result[0]):
                    matrix[result[1][i]][0] += d

            else:
                raise NotImplementedError(
                    f'Class {self.__class__.__name__} does not implement {method_name}')

        return matrix

    @property
    def dimension(self):
        return self.mesh.elements_groups[0].elements[0].dimension

    def get_row_col_indices(self, element):

        indexes = [self.mesh.node_to_index[point] for point in element.points]

        # positions = finite_elements.core.global_matrix_positions(dimension=self.dimension,
        #                                                           nodes_number=len(self.mesh.nodes))
        positions = self.positions

        row_ind, col = [], []
        for index in indexes:
            for i in range(element.dimension):
                row_ind.extend(len(indexes) * element.dimension * [positions[(index, i + 1)]])
                col.append(positions[(index, i + 1)])

        # col_ind = []
        # for index in indexes:
        #     for i in range(element.dimension):
        #         col_ind.extend(col)

        col_ind = col * (len(indexes) * element.dimension)

        return row_ind, col_ind

    def get_source_matrix_length(self):
        return len(self.mesh.nodes) * self.dimension + len(self.continuity_conditions) \
            + len(self._boundary_conditions)

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

        permeabilities = self.elements_permeability()
        mu_max = max(permeabilities)
        mu_min = min(permeabilities)
        colors = []
        for elements_group in self.mesh.elements_groups:
            x = (elements_group.elements[0].mu_total - mu_min) / (mu_max - mu_min)
            color = (color_map[0][0] - (color_map[0][0] - color_map[1][0]) * x,
                     color_map[0][1] - (color_map[0][1] - color_map[1][1]) * x,
                     color_map[0][2] - (color_map[0][2] - color_map[1][2]) * x)
            colors.append(color)

        norm = mpl.colors.Normalize(vmin=mu_min, vmax=mu_max)
        sm = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)  # plt.get_cmap('jet_r')
        sm.set_array([])
        cbar = fig.colorbar(sm, ticks=npy.linspace(mu_min, mu_max, 10))
        cbar.set_label('permeability')
        for i, elements_group in enumerate(self.mesh.elements_groups):
            for element in elements_group.elements:
                element.plot(ax=ax, color=colors[i], fill=True)

        return ax

    def plot_continuity_condition(self, ax=None):
        if ax is None:
            ax = self.mesh.plot()

        for i, continuity_condition in enumerate(self.continuity_conditions):
            continuity_condition.node1.MPLPlot(ax=ax, color='C{}'.format(i % 10))
            continuity_condition.node2.MPLPlot(ax=ax, color='C{}'.format(i % 10))
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
