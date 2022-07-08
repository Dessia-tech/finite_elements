#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Library: finite_elements (core.py)
"""

# import matplotlib as mpl
# import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as npy
import matplotlib.tri as mtri
# import volmdlr as vm
# import volmdlr.mesh as vmmesh
import math
# from scipy import sparse
# from scipy import linalg
# import time 
# from dessia_common import DessiaObject
# from typing import TypeVar, List, Tuple


cdict = {'red':  [(0.0, 0.0, 0.0),
                   (1.0, 1.0, 1.0)],
         'green': [(0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)],
         'blue':  [(0.0, 1.0, 1.0),
                   (1.0, 0.0, 0.0)]}

blue_red = LinearSegmentedColormap('BLueRed', cdict)

MU = 4*math.pi*1e-7


def get_bmin_bmax(Bs, Bmax=None, Bmin=None):
    if Bmax is None and Bmin is None:
        B_max, B_min = max(Bs), min(Bs)
    elif Bmax is not None and Bmin is None:
        B_max, B_min = Bmax, min(Bs)
    elif Bmax is None and Bmin is not None:
        B_max, B_min = max(Bs), Bmin
    else:
        B_max, B_min = Bmax, Bmin

    return B_max, B_min

def get_colors(Bs, B_max=None, B_min=None):
    color_map = ((0,0,1), (1,0,0))

    B_to_color = {}
    for B in Bs:
        if B > B_max:
            x = 1
        else:
            x = (B - B_min) / (B_max - B_min)

        color = (color_map[0][0]-(color_map[0][0]-color_map[1][0])*x,
                 color_map[0][1]-(color_map[0][1]-color_map[1][1])*x,
                 color_map[0][2]-(color_map[0][2]-color_map[1][2])*x)
        B_to_color[B] = color

    return B_to_color

def global_matrix_positions(dimension, nodes_number):

    positions = {}
    count = 0
    for i in range(nodes_number):
        for j in range(dimension):
            positions[(i, j+1)] = count
            count += 1

    return positions

def get_triangulation(mesh):
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

    x = npy.asarray(x_list)
    y = npy.asarray(y_list)
    triang = mtri.Triangulation(x, y, triangles)

    return triang
