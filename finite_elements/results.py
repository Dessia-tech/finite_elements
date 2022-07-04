#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing objects related to finite elements analysis results
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
# from matplotlib.colors import LinearSegmentedColormap
import numpy as npy
# import matplotlib.tri as mtri
import volmdlr as vm
import volmdlr.mesh as vmmesh
# import math
# from scipy import sparse
# from scipy import linalg
# import time 
from dessia_common import DessiaObject
from typing import List #Tuple, TypeVar
from finite_elements.core import MU, blue_red
import finite_elements.core
from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner
import matplotlib.cm as cm


class Result(DessiaObject):
    """
    :param mesh: The mesh on which the result has been solved.
    :type mesh: a Mesh object
    :param result_vector: The solution vector for the magnetic potential A.
    :type result_vector: a list of float
    """
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True
    def __init__(self, mesh: vmmesh.Mesh, result_vector: List[float]):
        self.mesh = mesh
        self.result_vector = result_vector
        
        DessiaObject.__init__(self, name='')
        
    def magnetic_field_per_element(self):
        element_to_magnetic_field = {}
        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:
                element_form_functions = element.form_functions
                indexes = [self.mesh.node_to_index[element.points[0]],
                           self.mesh.node_to_index[element.points[1]],
                           self.mesh.node_to_index[element.points[2]]]
                b1 = element_form_functions[0][1]
                c1 = element_form_functions[0][2]
                b2 = element_form_functions[1][1]
                c2 = element_form_functions[1][2]
                b3 = element_form_functions[2][1]
                c3 = element_form_functions[2][2]
                B_x = float(c1*self.result_vector[indexes[0]] + c2*self.result_vector[indexes[1]] + c3*self.result_vector[indexes[2]])
                B_y = float(-b1*self.result_vector[indexes[0]] - b2*self.result_vector[indexes[1]] - b3*self.result_vector[indexes[2]])
                element_to_magnetic_field[element] = vm.Vector2D(B_x, B_y)
        return element_to_magnetic_field

    def magnetic_field_norm(self):
        element_to_magnetic_field = self.magnetic_field_per_element()
        Bs = [B.norm() for B in list(element_to_magnetic_field.values())]

        return Bs

    def torque(self, air_gap_elements_group_name, length_motor, radius_stator, radius_rotor, nb_notches):
        """ 
        Computes the resistant magnetic torque when the rotor is blocked and \
        the current inside the stator is evolving. Unit : N.m.
        Uses Arkkio's method, based on the Maxwell Stress Tensor.
        
        :param air_gap_elements_group_name: The name given to the gap's ElementsGroup.
        :type air_gap_elements_group_name: str
        :param length_motor: The length of the Machine.
        :type length_motor: float
        :param radius_stator: The inner radius of the Stator.
        :type radius_stator: float
        :param radius_rotor: The outter radius of the Rotor.
        :type radius_rotor: float
        :param nb_notches: The number of notches of the Stator.
        :type nb_notches: int
        """
        element_to_magnetic_field = self.magnetic_field_per_element()
        
        for elements_group in self.mesh.elements_groups:
            if elements_group.name == air_gap_elements_group_name:
                gap_elements_group = elements_group
                break
        
        somme = 0
        i = 0
        
#        r = (radius_stator - radius_rotor)/2 + radius_rotor
        
#        fig, ax = plt.subplots()
        
        for element in gap_elements_group.elements:
            vector_B = element_to_magnetic_field[element]
            
            element_center = element.center
            e_r = vm.Vector2D(element_center.vector)
            e_r.Normalize()
            e_teta = vm.Vector2D((-e_r[1], e_r[0]))
            B_r = vector_B.Dot(e_r)
            B_teta = vector_B.Dot(e_teta)
            r_Br_Bteta = element_center.Norm() * B_r * B_teta
#            r_Br_Bteta = r * B_r * B_teta
            dS = element.area 
            
#            e_r.plot(ax=ax, origin=element_center, amplitude=0.005, color='b')
#            e_teta.plot(ax=ax, origin=element_center, amplitude=0.005, color='g')
#            vector_B.plot(ax=ax, origin=element_center, amplitude=0.005, color='r')
            
            somme += r_Br_Bteta * dS
            i += 1
        
#        print('nb elements in airgap', i)
        T = length_motor/(MU * (radius_stator - radius_rotor)) * somme
        
        return T 
    
    def maxwell_stress_tensor(self, element):
        """ 
        Computes the Maxwell stress tensor for one element.
        Returns a list made of sigma_r_r, sigma_r_theta and sigma_theta_theta, \
        since sigma_r_theta = sigma_theta_r.
        
        :param element: The element on which the tensor is computed.
        :type element: a TriangularElement object
        """
        element_to_magnetic_field = self.magnetic_field_per_element()
        vector_B = element_to_magnetic_field[element]
        element_center = element.center
        e_r = vm.Vector2D(element_center.vector)
        e_r.Normalize()
        e_teta = vm.Vector2D((-e_r[1], e_r[0]))
        B_r = vector_B.Dot(e_r)
        B_teta = vector_B.Dot(e_teta)
        
        sigma_rr = 1/MU * B_r**2 - 1/(2*MU) * vector_B.Norm()**2
        sigma_rteta = 1/MU * B_r * B_teta
        sigma_tetateta = 1/MU * B_teta**2 - 1/(2*MU) * vector_B.Norm()**2
        sigma_rr_rteta_tetateta = [sigma_rr, sigma_rteta, sigma_tetateta]
        return sigma_rr_rteta_tetateta

    def brbtetha(self, gap_elements_group):
        element_to_magnetic_field = self.magnetic_field_per_element()
        all_BrBtetha = []
        for element in gap_elements_group.elements:
            vector_B = element_to_magnetic_field[element]
            
            element_center = element.center
            e_r = vm.Vector2D(element_center.vector)
            e_r.Normalize()
            e_teta = vm.Vector2D((-e_r[1], e_r[0]))
            B_r = vector_B.Dot(e_r)
            B_teta = vector_B.Dot(e_teta)
            
            all_BrBtetha.append(B_r * B_teta)
        return all_BrBtetha

    def plot_brbtetha(self, ax=None, air_gap_elements_group_name='Gap ring'):
        if ax is None:
            fig, ax = plt.subplots()
            # ax = self.mesh.plot()
        else:
            fig = plt.gcf()
            # self.mesh.plot(ax=ax)

        for elements_group in self.mesh.elements_groups:
            if elements_group.name == air_gap_elements_group_name:
                gap_elements_group = elements_group
                break

        all_BrBtetha = self.brbtetha(gap_elements_group)

        # color_map = ((0,0,1), (1,0,0))
        jet = plt.get_cmap('jet')
        # Bs = [B.Norm() for B in list(element_to_magnetic_field.values())]
        
        BrBtetha_max, BrBtetha_min = max(all_BrBtetha), min(all_BrBtetha)
        
        B_to_color = {}
        all_colors = []
        for B in all_BrBtetha:
            if B > BrBtetha_max:
                x = 1
            else:
                x = (B - BrBtetha_min) / (BrBtetha_max - BrBtetha_min)
            # color = (color_map[0][0]-(color_map[0][0]-color_map[1][0])*x, 
            #          color_map[0][1]-(color_map[0][1]-color_map[1][1])*x,
            #          color_map[0][2]-(color_map[0][2]-color_map[1][2])*x)
            color = jet(int(x*256))[:3]
            
            B_to_color[B] = color
            all_colors.append(color)
        # print(B_to_color)
        for i, element in enumerate(gap_elements_group.elements):
            # element.plot(ax=ax, color=B_to_color[element_to_magnetic_field[element].Norm()], fill=True)
            element.plot(ax=ax, color=all_colors[i], fill=True)
            
        
        norm = mpl.colors.Normalize(vmin=BrBtetha_min, vmax=BrBtetha_max)
        sm = plt.cm.ScalarMappable(cmap=jet, norm=norm) 
        
        sm.set_array([])
        cbar = fig.colorbar(sm, ticks=npy.linspace(BrBtetha_min, BrBtetha_max, 10))
        # cbar = fig.colorbar(sm, ticks=npy.linspace(-0.9, 0.8, 10))
        cbar.set_label('Br*Btetha')
        
        return ax

    def plot_magnetic_field_vectors(self, ax=None, amplitude=0.005,
                                    Bmax=None, Bmin=None):
        """
        Plots the mesh with a field of vectors representing the induction field \
        and its direction inside the Machine. 
        """
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()

        element_to_magnetic_field = self.magnetic_field_per_element()
        Bs = self.magnetic_field_norm()

        B_max, B_min = finite_elements.core.get_bmin_bmax(Bs, Bmin, Bmax)
        B_to_color = finite_elements.core.get_colors(Bs, B_max=B_max, B_min=B_min)

        for element, B in element_to_magnetic_field.items():
            B.plot(amplitude=amplitude, origin=element.center, ax=ax,
                   color=B_to_color[B.norm()], normalize=True)

        norm = mpl.colors.Normalize(vmin=B_min, vmax=B_max)
        sm = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ticks=npy.linspace(B_min, B_max, 10))
        cbar.set_label('Magnetic Field in Tesla')

        return ax
    
    def plot_magnetic_field(self, ax=None, Bmax=None):
        """
        Plots the mesh with colored triangular elements representing the \
        intensity of the induction field inside the Machine. 
        """
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()
        element_to_magnetic_field = self.magnetic_field_per_element()
        Bs =  self.magnetic_field_norm()
        
        B_max, B_min = finite_elements.core.get_bmin_bmax(Bs, Bmin=None, Bmax=Bmax)
        B_to_color = finite_elements.core.get_colors(Bs, B_max=B_max, B_min=B_min)

        for group in self.mesh.elements_groups:
            for element in group.elements:
                element.plot(ax=ax, color=B_to_color[element_to_magnetic_field[element].norm()], fill=True)
        
        norm = mpl.colors.Normalize(vmin=B_min, vmax=B_max)
        sm = plt.cm.ScalarMappable(cmap=blue_red, norm=norm) 
        sm.set_array([])
        cbar = fig.colorbar(sm, ticks=npy.linspace(B_min, B_max, 10))
        cbar.set_label('Magnetic Field in Tesla')
        
        return ax
    
    def plot_magnetic_field_contour(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()
        ax.set_aspect('equal')

        element_to_magnetic_field = self.magnetic_field_per_element()

        x = []
        y = []
        Z = []
        for group in self.mesh.elements_groups:
            for element in group.elements:
                x_center, y_center = element.center
                x.append(x_center)
                y.append(y_center)
                Z.append(element_to_magnetic_field[element].norm())

        tri = Triangulation(x, y)

        #-----------------------------------------------------------------------------
        # Improving the triangulation before high-res plots: removing flat triangles
        #-----------------------------------------------------------------------------
        # masking badly shaped triangles at the border of the triangular mesh.
        min_circle_ratio = -1 
        mask = TriAnalyzer(tri).get_flat_tri_mask(min_circle_ratio)
        tri.set_mask(mask)

        # refining the data
        refiner = UniformTriRefiner(tri)
        subdiv = 3
        tri_refi, z_test_refi = refiner.refine_field(Z, subdiv=subdiv)

        levels = npy.arange(0., 1., 0.05)
        cmap = cm.get_cmap(name='Blues', lut=None)
        ax.tricontour(tri_refi, z_test_refi, levels=levels, #cmap=cmap,
                      linewidths=[2.0, 0.5, 1.0, 0.5])
#        ax.triplot(tri_refi, color='0.97')
#        ax.triplot(tri, color='0.7')
#        ax.tricontour(x, y, Z)
        return ax 

    def plot_potential_vector(self, ax=None):
        """
        Plots the mesh with colored triangular elements representing the \
        intensity of the magnetic potential inside the Machine. 
        """

        triang = finite_elements.core.get_triangulation(self.mesh)
        z = npy.asarray([p for p in self.result_vector[:len(self.mesh.nodes)]]) #p[0]
        z_min, z_max = min(z), max(z)

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()
        ax.tricontourf(triang, z, cmap=blue_red)
        ax.triplot(triang, 'k-')
        ax.set_title('Triangular grid')

        norm = mpl.colors.Normalize(vmin=z_min,vmax=z_max)
        sm = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ticks=npy.linspace(z_min, z_max, 10))
        cbar.set_label('Potential Vector')

        return ax