#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing objects related to different conditions types
"""

import volmdlr as vm
import volmdlr.mesh as vmmesh
# import math
# from scipy import sparse
# from scipy import linalg
# import time 
from dessia_common import DessiaObject
# from typing import List


class BoundaryCondition(DessiaObject):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, application, value: float, dimension, name: str = ''):
        self.application = application
        self.value = value
        self.dimension = dimension

        DessiaObject.__init__(self, name='')


class NodeBoundaryCondition(BoundaryCondition):
    # _standalone_in_db = False
    # _non_serializable_attributes = []
    # _non_eq_attributes = ['name']
    # _non_hash_attributes = ['name']
    # _generic_eq = True

    def __init__(self, application, value: float, dimension, name: str = ''):
        self.application = application #Node
        self.value = value
        self.dimension = dimension

        BoundaryCondition.__init__(self, application, value, dimension, name='')

class ElementBoundaryCondition(BoundaryCondition):
    # _standalone_in_db = False
    # _non_serializable_attributes = []
    # _non_eq_attributes = ['name']
    # _non_hash_attributes = ['name']
    # _generic_eq = True

    def __init__(self, application, value: float, dimension, name: str = ''):
        self.application = application #Element
        self.value = value
        self.dimension = dimension

        BoundaryCondition.__init__(self, application, value, dimension, name='')


class ContinuityCondition(DessiaObject):
    """ 
    The continuity conditions link the value of vector potential A between two \
    nodes. It is used to describe periodic or antiperiodic conditions inside the \
    mesh. 
    
    :param node1: The first node.
    :type node1: volmdlr.Point2D object
    :param node2: The second node.
    :type node2: volmdlr.Point2D object
    :param value: If value is equal to 1, the continuity condition is periodic. \
    If value is equal to -1, the continuity condition is antiperiodic. A(node1) = value*A(node2).
    :type value: int
    """
    def __init__(self, node1: vm.Point2D, node2: vm.Point2D, value):
        self.node1 = node1
        self.node2 = node2 
        self.value = value
        
        DessiaObject.__init__(self, name='')