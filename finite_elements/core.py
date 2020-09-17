#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 11:58:50 2020

@author: gasmi
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as npy
import matplotlib.tri as mtri
import volmdlr as vm
import volmdlr.mesh as vmmesh
import math
from scipy import sparse
from scipy import linalg
import time  
from dessia_common import DessiaObject
from typing import TypeVar, List, Tuple,Dict

steel=[210,0.3,7800*1E3]
aluminium=[69,0.346]

class Materials(DessiaObject):
    """
    Sets the physicals constants of the material (of the element group)
    :param materials: Dictionnary of elements groups and their respective young module and poisson coefficient 
    :type materials:  Dictionnary of elements groups and a list of their young modules and poisson coefficients
    """
    def __init__(self,materials_properties:Dict[vmmesh.ElementsGroup,List[float]]):
        self.materials_properties=materials_properties
        
