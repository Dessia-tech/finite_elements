#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing objects related to finite elements analysis
"""

from typing import List
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as npy
import volmdlr as vm
import volmdlr.mesh as vmmesh
# import time
import scipy
from scipy import sparse
import scipy.sparse.linalg
from scipy.sparse import csc_matrix

from dessia_common.core import DessiaObject

import finite_elements.elements
from finite_elements.loads import ElementsLoad, EdgeLoad, NodeLoad, MagnetLoad
from finite_elements.conditions import ContinuityCondition
from finite_elements.results import Result
from finite_elements.core import blue_red


def node_boundary_conditions_to_dict(node_boundary_conditions):
    """
    Creates

    :param node_boundary_conditions: DESCRIPTION
    :type node_boundary_conditions: TYPE

    :return: DESCRIPTION
    :rtype: TYPE
    """

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
    """
    Creates

    :param node_boundary_conditions_dict: DESCRIPTION
    :type node_boundary_conditions_dict: TYPE

    :return: DESCRIPTION
    :rtype: TYPE
    """

    node_boundary_conditions = []
    for key, value in node_boundary_conditions_dict.items():
        node_boundary_conditions.append(
            finite_elements.conditions.NodeBoundaryCondition(
                application=key[0],
                value=value,
                dimension=key[1]))
    return node_boundary_conditions


def node_loads_to_dict(node_loads):
    """
    Creates

    :param node_loads: DESCRIPTION
    :type node_loads: TYPE

    :return: DESCRIPTION
    :rtype: TYPE
    """

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
    """
    Creates

    :param node_loads_dict: DESCRIPTION
    :type node_loads_dict: TYPE

    :return: DESCRIPTION
    :rtype: TYPE
    """

    node_loads = []
    for key, value in node_loads_dict.items():
        node_loads.append(NodeLoad(
            node=key[0],
            value=value,
            dimension=key[1]))
    return node_loads


class FiniteElements(DessiaObject):
    """
    This class

    :param mesh: DESCRIPTION
    :type mesh: vmmesh.Mesh
    :param element_loads: DESCRIPTION
    :type element_loads: List[ElementsLoad]
    :param edge_loads: DESCRIPTION
    :type edge_loads: List[EdgeLoad]
    :param node_loads: DESCRIPTION
    :type node_loads: List[NodeLoad]
    :param magnet_loads: DESCRIPTION
    :type magnet_loads: List[MagnetLoad]
    :param continuity_conditions: DESCRIPTION
    :type continuity_conditions: List[ContinuityCondition]
    :param node_boundary_conditions: DESCRIPTION
    :type node_boundary_conditions: List[finite_elements.conditions.NodeBoundaryCondition]
    :param edge_boundary_conditions: DESCRIPTION
    :type edge_boundary_conditions: List[finite_elements.conditions.EdgeBoundaryCondition]
    :param element_boundary_conditions: DESCRIPTION
    :type element_boundary_conditions: List[finite_elements.conditions.ElementBoundaryCondition]
    :param plane_strain: DESCRIPTION, defaults to None
    :type plane_strain: bool, optional
    :param plane_stress: DESCRIPTION, defaults to None
    :type plane_stress: bool, optional
    """

    def __init__(self, mesh: vmmesh.Mesh,
                 element_loads: List[ElementsLoad],
                 edge_loads: List[EdgeLoad],
                 node_loads: List[NodeLoad],
                 magnet_loads: List[MagnetLoad],
                 continuity_conditions:
                     List[ContinuityCondition],
                 node_boundary_conditions:
                     List[finite_elements.conditions.NodeBoundaryCondition],
                 edge_boundary_conditions:
                     List[finite_elements.conditions.EdgeBoundaryCondition],
                 element_boundary_conditions:
                     List[finite_elements.conditions.ElementBoundaryCondition],
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
        """
        Creates

        :return: DESCRIPTION
        :rtype: TYPE
        """

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
        """
        Creates

        :return: DESCRIPTION
        :rtype: TYPE
        """

        element_to_node_boundary_conditions = []
        for element_condition in self.element_boundary_conditions:

            matrix_factors = element_condition.application.element_to_node_factors()

            for p_index, point in enumerate(element_condition.application.points):
                element_to_node_boundary_conditions.append(
                    finite_elements.conditions.NodeBoundaryCondition(
                        application=point,
                        value=element_condition.value * matrix_factors[p_index],
                        dimension=element_condition.dimension))
        return element_to_node_boundary_conditions

    def boundary_conditions_edge_to_node(self):
        """
        Creates

        :return: DESCRIPTION
        :rtype: TYPE
        """

        edge_to_node_boundary_conditions = []
        for edge_condition in self.edge_boundary_conditions:
            for point in [edge_condition.application.start,
                          edge_condition.application.end]:
                edge_to_node_boundary_conditions.append(
                    finite_elements.conditions.NodeBoundaryCondition(
                        application=point,
                        value=edge_condition.value * 0.5,
                        dimension=edge_condition.dimension))
        return edge_to_node_boundary_conditions

    def c_matrix_boundary_conditions(self):
        """
        Creates

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if self._boundary_conditions:
            node_boundary_conditions = self._boundary_conditions
        else:
            node_boundary_conditions = self.node_boundary_conditions[:]
            # element to node
            node_boundary_conditions.extend(self.boundary_conditions_element_to_node())
            # edge to node
            node_boundary_conditions.extend(self.boundary_conditions_edge_to_node())

            # node_bc to_dict
            node_boundary_conditions_dict = node_boundary_conditions_to_dict(
                node_boundary_conditions)

            # node_bc dict from_dict
            node_boundary_conditions = node_boundary_from_dict(
                node_boundary_conditions_dict)
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
    def dimension(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return self.mesh.elements_groups[0].elements[0].dimension

    @property
    def elements_name(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return self.mesh.elements_groups[0].elements[0].__class__.__name__

    def elements_permeability(self):
        """
        Creates

        :return: DESCRIPTION
        :rtype: TYPE
        """

        permeabilities = [elements_group.elements[0].mu_total for elements_group in
                          self.mesh.elements_groups]

        return permeabilities

    def k_matrix(self, method_name):
        if method_name == 'dense':
            return self.k_matrix_dense()
        if method_name == 'sparse':
            return self.k_matrix_sparse()
        raise NotImplementedError(
            f'Class {self.__class__.__name__} does not implement {method_name} k matrix')

    def k_matrix_data(self):
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

    def k_matrix_dense(self):

        return self.matrix_dense(method_name='k_matrix_data')

    def k_matrix_sparse(self):

        return self.matrix_sparse(method_name='k_matrix_data')

    def m_matrix(self, method_name):
        if method_name == 'dense':
            return self.m_matrix_dense()
        if method_name == 'sparse':
            return self.m_matrix_sparse()
        raise NotImplementedError(
            f'Class {self.__class__.__name__} does not implement {method_name} m matrix')

    def m_matrix_data(self):
        row_ind, col_ind, data = [], [], []
        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:
                data.extend(element.elementary_mass_matrix())
                row_ind_n, col_ind_n = self.get_row_col_indices(element)
                row_ind.extend(row_ind_n)
                col_ind.extend(col_ind_n)
        return data, row_ind, col_ind

    def m_matrix_dense(self):

        return self.matrix_dense(method_name='m_matrix_data')

    def m_matrix_sparse(self):

        return self.matrix_sparse(method_name='m_matrix_data')

    def matrix_dense(self, method_name):
        matrix = npy.zeros((len(self.mesh.nodes) * self.dimension,
                            len(self.mesh.nodes) * self.dimension))
        if hasattr(self, method_name):
            data, row_ind, col_ind = getattr(self, method_name)()
            for i, data_i in enumerate(data):
                matrix[row_ind[i], col_ind[i]] += data_i
        else:
            raise NotImplementedError(
                f'Class {self.__class__.__name__} does not implement {method_name}')
        return matrix

    def matrix_sparse(self, method_name):
        if hasattr(self, method_name):
            data, row_ind, col_ind = getattr(self, method_name)()
            dict_data = {}
            for i, data_i in enumerate(data):
                try:
                    dict_data[(row_ind[i], col_ind[i])] += data_i
                except KeyError:
                    dict_data[(row_ind[i], col_ind[i])] = data_i
            data_new, row_ind_new, col_ind_new = [], [], []
            for key, value in dict_data.items():
                data_new.append(value)
                row_ind_new.append(key[0])
                col_ind_new.append(key[1])
            matrix = sparse.csc_matrix((data_new, (row_ind_new, col_ind_new)))
        else:
            raise NotImplementedError(
                f'Class {self.__class__.__name__} does not implement {method_name}')
        return matrix

    def loads_element_to_node(self):
        """
        Creates

        :return: DESCRIPTION
        :rtype: TYPE
        """

        element_to_node_loads = []
        for elements_load in self.element_loads:
            loads_per_element = elements_load.value_per_element

            for j, element in enumerate(elements_load.elements):

                matrix_factors = element.element_to_node_factors()

                for p_index, point in enumerate(element.points):
                    element_to_node_loads.append(NodeLoad(
                        node=point,
                        value=loads_per_element[j] * matrix_factors[p_index],
                        dimension=elements_load.dimension))

        return element_to_node_loads

    def loads_edge_to_node(self):
        """
        Creates

        :return: DESCRIPTION
        :rtype: TYPE
        """

        edge_to_node_loads = []
        for edge_load in self.edge_loads:
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
        """
        Creates

        :return: DESCRIPTION
        :rtype: TYPE
        """

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
        for _, load in enumerate(node_loads):
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
        """
        Creates

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if self._boundary_conditions:
            node_boundary_conditions = self._boundary_conditions
        else:
            node_boundary_conditions = self.node_boundary_conditions[:]
            # element to node
            node_boundary_conditions.extend(self.boundary_conditions_element_to_node())
            # edge to node
            node_boundary_conditions.extend(self.boundary_conditions_edge_to_node())

            # node_bc to_dict
            node_boundary_conditions_dict = node_boundary_conditions_to_dict(
                node_boundary_conditions)

            # node_bc dict from_dict
            node_boundary_conditions = node_boundary_from_dict(
                node_boundary_conditions_dict)

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
        """
        Creates

        :return: DESCRIPTION
        :rtype: TYPE
        """

        data, row_ind = [], []
        for magnet_load in self.magnet_loads:
            for linear_element in magnet_load.contour_linear_elements():
                indexes = [self.mesh.node_to_index[linear_element.points[0]],
                           self.mesh.node_to_index[linear_element.points[1]]]
                length = linear_element.length()
                dl_parameter = vm.Vector2D([-linear_element.interior_normal[1],
                                            linear_element.interior_normal[0]])
                data.append(magnet_load.magnetization_vector.Dot(dl_parameter) * length / 2)
                data.append(magnet_load.magnetization_vector.Dot(dl_parameter) * length / 2)
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
    #               element_loads: List[ElementsLoad],
    #               node_loads: List[NodeLoad],
    #               magnet_loads: List[MagnetLoad],
    #               continuity_conditions: List[ContinuityCondition],
    #               node_boundary_conditions:
    #                   List[finite_elements.conditions.NodeBoundaryCondition],
    #               element_boundary_conditions:
    #                   List[finite_elements.conditions.ElementBoundaryCondition]):
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
        """
        Creates

        :return: DESCRIPTION
        :rtype: TYPE
        """

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

        method_names = ['k_matrix_data', 'c_matrix_continuity_conditions',
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
        """
        Creates

        :return: DESCRIPTION
        :rtype: TYPE
        """

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
                for i, data in enumerate(result[0]):
                    matrix[result[1][i]][0] += data

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

    def modal_analysis(self, order, k):
        if order == 'largest':
            # REMARKS
            # eigh & eigsh had been compared => eigsh is the best
            # dim=4850, k=30: eigsh: 3.973 | eigh => 13.139
            # dim=14949, k=30: eigsh: 2.449 | eigh => 293.837

            # =============================================================================
            # ### dense
            # =============================================================================
            # dimension = len(self.mesh.nodes)*self.dimension
            # t = time.time()
            # m_matrix_dense = self.m_matrix_dense()
            # k_matrix_dense = self.k_matrix_dense()
            # print('dense')
            # print('k, m matrices => ', time.time()-t)
            # t = time.time()
            # eigvals, eigvecs = scipy.linalg.eigh(k_matrix_dense,
            #                                     m_matrix_dense,
            #                                     check_finite = False,
            #                                     eigvals = (dimension - k,
            #                                                dimension - 1))
            # print('eigh => ', time.time()-t)
            # print('************************')

            # =============================================================================
            # ### sparse
            # =============================================================================
            # REMARKS
            # self.m_matrix_sparse() is faster than csc_matrix(self.m_matrix_dense())

            # t = time.time()
            m_matrix_sparse = self.m_matrix_sparse()
            k_matrix_sparse = self.k_matrix_sparse()
            # print('sparse')
            # print('k, m matrices => ', time.time()-t)

            # t = time.time()
            eigvals, eigvecs = scipy.sparse.linalg.eigsh(A=k_matrix_sparse,
                                                         M=m_matrix_sparse,
                                                         which='LM',
                                                         k=k)
            # print('eigsh => ', time.time()-t)
            # print('************************')

            return eigvals, eigvecs.T

        if order == 'smallest':
            m_matrix_sparse = self.m_matrix_sparse()
            k_matrix_sparse = csc_matrix(scipy.linalg.inv(self.k_matrix_dense()))
            eigvals, eigvecs = scipy.sparse.linalg.eigs(k_matrix_sparse @ m_matrix_sparse,
                                                        which='LM',
                                                        k=k)
            return 1 / eigvals, eigvecs.T

        raise ValueError("Order parameter should be either 'largest' or 'smallest'")

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

        k_sparse = self.create_matrix()
        force = self.create_source_matrix()

        try:
            x_solution = sparse.linalg.spsolve(k_sparse, force,
                                               permc_spec='NATURAL',
                                               use_umfpack=True)

        except sparse.linalg.MatrixRankWarning as exc:
            print('MatricRankWarning')
            raise NotImplementedError from exc

        x_solution = list(x_solution)
        # print('apres')
        return Result(self.mesh, x_solution)

    def plot_elements_loads(self, ax=None):
        """
        Plots the mesh. The triangular elements are filled red if they are a
            source of magnetic potential thanks to their current density

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if ax is None:
            _, ax = plt.subplots()
        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:
                element.plot(ax=ax,
                             edge_style=vm.core.EdgeStyle(color='w'),
                             fill_color='w',
                             fill=True)

        for load in self.element_loads:
            for element in load.elements:
                element.plot(ax=ax,
                             edge_style=vm.core.EdgeStyle(color='r'),
                             fill_color='r',
                             fill=True)
        return ax

    def plot_elements_permeability(self, ax=None):
        """
        Plots the mesh with colored triangular elements depending on the
            permeability of the material the element meshes

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
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
            x_param = (elements_group.elements[0].mu_total - mu_min) / (mu_max - mu_min)
            color = (color_map[0][0] - (color_map[0][0] - color_map[1][0]) * x_param,
                     color_map[0][1] - (color_map[0][1] - color_map[1][1]) * x_param,
                     color_map[0][2] - (color_map[0][2] - color_map[1][2]) * x_param)
            colors.append(color)

        norm = mpl.colors.Normalize(vmin=mu_min, vmax=mu_max)
        sm_param = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)  # plt.get_cmap('jet_r')
        sm_param.set_array([])
        cbar = fig.colorbar(sm_param, ticks=npy.linspace(mu_min, mu_max, 10))
        cbar.set_label('permeability')
        for i, elements_group in enumerate(self.mesh.elements_groups):
            for element in elements_group.elements:
                element.plot(ax=ax,
                             edge_style=vm.core.EdgeStyle(color=colors[i]),
                             fill_color=colors[i],
                             fill=True)

        return ax

    def plot_continuity_condition(self, ax=None):
        """
        Plots

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if ax is None:
            ax = self.mesh.plot()

        for i, continuity_condition in enumerate(self.continuity_conditions):
            continuity_condition.node1.MPLPlot(ax=ax, color=f'C{i % 10}')
            continuity_condition.node2.MPLPlot(ax=ax, color=f'C{i % 10}')
        return ax

    def plot_magnet_loads(self, ax=None):
        """
        Plots the mesh. The triangular elements are filled red if they are a
            source of magnetic potential thanks to their magnetization. The contour
            of the magnetizating mesh is drawn in blue.

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if ax is None:
            _, ax = plt.subplots()
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
