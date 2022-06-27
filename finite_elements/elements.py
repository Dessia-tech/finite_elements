#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing objects related to different finite elements types
"""

# import matplotlib as mpl
# import matplotlib.pyplot as plt
# from matplotlib.colors import LinearSegmentedColormap
# import numpy as npy
# import matplotlib.tri as mtri
# import volmdlr as vm
import volmdlr.mesh as vmmesh
# import math
# from scipy import sparse
# from scipy import linalg
# import time 
# from dessia_common import DessiaObject
from typing import List #Tuple, TypeVar


class MagneticElement2D(vmmesh.TriangularElement2D):
    # _standalone_in_db = False
    # _non_serializable_attributes = []
    # _non_eq_attributes = ['name']
    # _non_hash_attributes = ['name']
    # _generic_eq = True
    def __init__(self, triangular_element: vmmesh.TriangularElement2D,
                 mu_total: float, name : str = ''):
        self.triangular_element = triangular_element
        vmmesh.TriangularElement2D.__init__(self, points=triangular_element.points)
        self.mu_total = mu_total

        # DessiaObject.__init__(self, name=name)

    def elementary_matrix(self):
        """
        Create the elementary matrix of the MagneticElement2D

        :return: (data, row_ind, col_ind)
        """

        element_form_functions = self.triangular_element.form_functions
        indexes = [self.mesh.node_to_index[self.triangular_element.points[0]],
                   self.mesh.node_to_index[self.triangular_element.points[1]],
                   self.mesh.node_to_index[self.triangular_element.points[2]]]
        b1 = element_form_functions[0][1]
        c1 = element_form_functions[0][2]
        b2 = element_form_functions[1][1]
        c2 = element_form_functions[1][2]
        b3 = element_form_functions[2][1]
        c3 = element_form_functions[2][2]

        row_ind = (indexes[0], indexes[0], indexes[0], indexes[1], indexes[1], indexes[1], indexes[2], indexes[2], indexes[2])
        col_ind = (indexes[0], indexes[1], indexes[2], indexes[0], indexes[1], indexes[2], indexes[0], indexes[1], indexes[2])
        data = (1/self.mu_total * (b1**2 + c1**2) * self.triangular_element.area,
                1/self.mu_total * (b1*b2 + c1*c2) * self.triangular_element.area,
                1/self.mu_total * (b1*b3 + c1*c3) * self.triangular_element.area,
                1/self.mu_total * (b1*b2 + c1*c2) * self.triangular_element.area,
                1/self.mu_total * (b2**2 + c2**2) * self.triangular_element.area,
                1/self.mu_total * (b2*b3 + c2*c3) * self.triangular_element.area,
                1/self.mu_total * (b1*b3 + c1*c3) * self.triangular_element.area,
                1/self.mu_total * (b2*b3 + c2*c3) * self.triangular_element.area,
                1/self.mu_total * (b3**2 + c3**2) * self.triangular_element.area)

        return (data, row_ind, col_ind)


class MagneticElementsGroup(vmmesh.ElementsGroup):
    # _standalone_in_db = False
    # _non_serializable_attributes = []
    # _non_eq_attributes = ['name']
    # _non_hash_attributes = ['name']
    # _generic_eq = True
    def __init__(self, magnetic_elements: List[MagneticElement2D],
                 mu_total: float, name: str):
        self.magnetic_elements = magnetic_elements
        self.triangular_elements = self._triangular_elements()
        vmmesh.ElementsGroup.__init__(self, elements=self.triangular_elements, name=name)
        self.mu_total = mu_total

        # DessiaObject.__init__(self, name=name)

    def _triangular_elements(self):
        return [element.triangular_element for element in self.magnetic_elements]

