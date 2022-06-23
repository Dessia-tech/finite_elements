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


class MagneticElementsGroup(vmmesh.ElementsGroup):
    # _standalone_in_db = False
    # _non_serializable_attributes = []
    # _non_eq_attributes = ['name']
    # _non_hash_attributes = ['name']
    # _generic_eq = True
    def __init__(self, elements: List[vmmesh.TriangularElement],
                 mu_total: float, name: str):
        vmmesh.ElementsGroup.__init__(self, elements=elements, name=name)
        self.mu_total = mu_total
        
        # DessiaObject.__init__(self, name=name)
