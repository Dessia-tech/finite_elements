#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Library: finite_elements (core.py)
"""

# import matplotlib as mpl
# import matplotlib.pyplot as plt
import math
from matplotlib.colors import LinearSegmentedColormap
import numpy as npy
import matplotlib.tri as mtri
# import volmdlr
# import volmdlr.mesh as vmmesh
# from scipy import sparse
# from scipy import linalg
# import time
from dessia_common import DessiaObject
# from typing import TypeVar, List, Tuple


cdict = {'red': [(0.0, 0.0, 0.0),
                 (1.0, 1.0, 1.0)],
         'green': [(0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)],
         'blue': [(0.0, 1.0, 1.0),
                  (1.0, 0.0, 0.0)]}

blue_red = LinearSegmentedColormap('BLueRed', cdict)

MU = 4 * math.pi * 1e-7


def get_bmin_bmax(param_bs, param_bmax=None, param_bmin=None):
    """

    :param param_bs: DESCRIPTION
    :type param_bs: TYPE
    :param param_bmax: DESCRIPTION, defaults to None
    :type param_bmax: TYPE, optional
    :param Bmin: DESCRIPTION, defaults to None
    :type Bmin: TYPE, optional

    :return: DESCRIPTION
    :rtype: TYPE
    """

    if param_bmax is None and param_bmin is None:
        param_bmax, param_bmin = max(param_bs), min(param_bs)
    elif param_bmax is not None and param_bmin is None:
        param_bmin = min(param_bs)
    elif param_bmax is None and param_bmin is not None:
        param_bmax = max(param_bs)

    return param_bmax, param_bmin


def get_colors(param_bs, param_bmax=None, param_bmin=None):
    """

    :param param_bs: DESCRIPTION
    :type param_bs: TYPE
    :param param_bmax: DESCRIPTION, defaults to None
    :type param_bmax: TYPE, optional
    :param param_bmin: DESCRIPTION, defaults to None
    :type param_bmin: TYPE, optional

    :return: DESCRIPTION
    :rtype: TYPE
    """

    color_map = ((0, 0, 1), (1, 0, 0))

    b_to_color = {}
    for param_b in param_bs:
        if param_b > param_bmax:
            x_param = 1
        else:
            x_param = (param_b - param_bmin) / (param_bmax - param_bmin)

        color = (color_map[0][0] - (color_map[0][0] - color_map[1][0]) * x_param,
                 color_map[0][1] - (color_map[0][1] - color_map[1][1]) * x_param,
                 color_map[0][2] - (color_map[0][2] - color_map[1][2]) * x_param)
        b_to_color[param_b] = color

    return b_to_color


def global_matrix_positions(dimension, nodes_number):
    """

    :param dimension: DESCRIPTION
    :type dimension: TYPE
    :param nodes_number: DESCRIPTION
    :type nodes_number: TYPE

    :return: DESCRIPTION
    :rtype: TYPE
    """

    positions = {}
    count = 0
    for i in range(nodes_number):
        for j in range(dimension):
            positions[(i, j + 1)] = count
            count += 1

    return positions


def get_triangulation(mesh):
    """

    :param mesh: DESCRIPTION
    :type mesh: TYPE

    :return: DESCRIPTION
    :rtype: TYPE
    """

    x_list, y_list = [], []
    for node in mesh.nodes:
        x_list.append(node[0])
        y_list.append(node[1])
    triangles = []
    for group in mesh.elements_groups:
        for element in group.elements:
            triangles.append([mesh.node_to_index[element.points[0]],
                              mesh.node_to_index[element.points[1]],
                              mesh.node_to_index[element.points[2]]])

    x_coor = npy.asarray(x_list)
    y_coor = npy.asarray(y_list)
    triang = mtri.Triangulation(x_coor, y_coor, triangles)

    return triang


class Material(DessiaObject):
    """
    This class defines a material with its physical characteristics

    :param elasticity_modulus: The unit of measurement of an object's resistance to being deformed
        elastically, when a stress is applied to it
    :type elasticity_modulus: float
    :param poisson_ratio:  The measure of the Poisson effect, the deformation (expansion or
       contraction) of a material in directions perpendicular to the specific direction of loading
    :type poisson_ratio: float
    :param mass_density: The substance's mass per unit of volume
    :type mass_density: float
    :param name: The name, defaults to ''
    :type name: str, optional
    """

    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self,
                 elasticity_modulus,
                 poisson_ratio,
                 mass_density,
                 name: str = ''):
        self.elasticity_modulus = elasticity_modulus
        self.poisson_ratio = poisson_ratio
        self.mass_density = mass_density
        self.name = name

        DessiaObject.__init__(self, name=name)
