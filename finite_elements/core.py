#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Library: finite_elements (core.py)
"""

# import matplotlib as mpl
# import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# import numpy as npy
# import matplotlib.tri as mtri
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
