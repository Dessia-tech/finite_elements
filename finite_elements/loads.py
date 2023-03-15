#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing objects related to different loads types
"""

from typing import List
import volmdlr as vm
import volmdlr.mesh as vmmesh
from dessia_common.core import DessiaObject


class ElementsLoad(DessiaObject):
    """
    Sets a load on the selected elements by imposing a source value for the
        current density vector J. Each element creates a source of the input value

    :param elements: The triangular elements.
    :type elements: List of volmdlr.TriangularElements objects
    :param value: Set the elements' current density vector J value.
    :type value: float
    """

    def __init__(self, elements: List[vmmesh.TriangularElement],
                 value: float,
                 dimension: int):
        self.elements = elements
        self.value = value
        self.dimension = dimension

        self.value_per_element = []
        total_area = sum(elem.area for elem in self.elements)
        for element in self.elements:
            self.value_per_element.append(value * element.area / total_area)
        DessiaObject.__init__(self, name='')


class ElementLoad(DessiaObject):
    """
    Sets a load on the selected elements by imposing a source value for the
        current density vector J. Each element creates a source of the input value

    :param elements: The triangular elements.
    :type elements: List of volmdlr.TriangularElements objects
    :param value: Set the elements' current density vector J value.
    :type value: float
    """

    def __init__(self, element: vmmesh.TriangularElement,
                 value: float,
                 dimension: int):
        self.element = element
        self.value = value
        self.dimension = dimension

        DessiaObject.__init__(self, name='')


class EdgeLoad(DessiaObject):
    """
    This class
    """

    def __init__(self, edge, value, dimension):
        self.edge = edge
        self.value = value
        self.dimension = dimension

        DessiaObject.__init__(self, name='')


class NodeLoad(DessiaObject):
    """
    Forces the value of the vector potential A at a node. To set a magnetic wall
        the value of the vector potential has to be set to A

    :param node: The node.
    :type node: volmdlr.Point2D object
    :param value: Set the node's vector potential A value.
    :type value: float
    """

    def __init__(self, node: vm.Point2D, value: float, dimension):
        self.node = node
        self.value = value
        self.dimension = dimension

        DessiaObject.__init__(self, name='')

    def c_matrix(self):
        """
        Defines
        """

        return ()

    def source_c_matrix(self):
        """
        Defines
        """

        return self.value


class MagnetLoad(DessiaObject):
    """
    Sets a load on the selected elements by imposing a source value for the
        magnetization vector M. Each element creates a source of the input value.

    :param elements: The triangular elements.
    :type elements: List of volmdlr.TriangularElements
    :param non_contour_nodes: Specify the nodes that are not part of the magnet's contour.
    :type non_contour_nodes: List of volmdlr.Point2D objects
    :param magnetization_vector: Set the elements' magnetization vector M.
    :type magnetization_vector: volmdlr.Vector2D object
    """

    def __init__(self, elements: List[vmmesh.TriangularElement],
                 non_contour_nodes: List[vm.Point2D],
                 magnetization_vector: vm.Vector2D):
        self.elements = elements
        self.non_contour_nodes = non_contour_nodes
        self.magnetization_vector = magnetization_vector

        self.element_magnetization_vector = magnetization_vector / len(elements)

        DessiaObject.__init__(self, name='')

    def contour_linear_elements(self):
        """
        Defines
        """

        linear_elements_count = {}
        for element in self.elements:
            for linear_element in element.linear_elements:
                if linear_element not in linear_elements_count:
                    linear_elements_count[linear_element] = 1
                else:
                    linear_elements_count[linear_element] += 1
        contour_linear_elements = []
        for linear_element, count in linear_elements_count.items():
            if count == 1 \
                and (linear_element.points[0] not in self.non_contour_nodes
                     or linear_element.points[1] not in self.non_contour_nodes):
                contour_linear_elements.append(linear_element)
        return contour_linear_elements
