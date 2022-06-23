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


class MagneticElement(vmmesh.TriangularElement2D):
    # _standalone_in_db = False
    # _non_serializable_attributes = []
    # _non_eq_attributes = ['name']
    # _non_hash_attributes = ['name']
    # _generic_eq = True
    def __init__(self, triangular_element: vmmesh.TriangularElement,
                 mu_total: float, name : str = ''):
        self.triangular_element = triangular_element
        vmmesh.TriangularElement2D.__init__(self, points=triangular_element.points)
        self.mu_total = mu_total

        # DessiaObject.__init__(self, name=name)

class MagneticElementsGroup(vmmesh.ElementsGroup):
    # _standalone_in_db = False
    # _non_serializable_attributes = []
    # _non_eq_attributes = ['name']
    # _non_hash_attributes = ['name']
    # _generic_eq = True
    def __init__(self, magnetic_elements: List[MagneticElement],
                 mu_total: float, name: str):
        self.magnetic_elements = magnetic_elements
        self.triangular_elements = self._triangular_elements()
        vmmesh.ElementsGroup.__init__(self, elements=self.triangular_elements, name=name)
        self.mu_total = mu_total

        # DessiaObject.__init__(self, name=name)

    def _triangular_elements(self):
        return [element.triangular_element for element in self.magnetic_elements]

