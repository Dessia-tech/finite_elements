#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Library: finite_elements (core.py)
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
from typing import TypeVar, List, Tuple

cdict = {'red':  [(0.0, 0.0, 0.0),
                   (1.0, 1.0, 1.0)],
         'green': [(0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)],
         'blue':  [(0.0, 1.0, 1.0),
                   (1.0, 0.0, 0.0)]}
blue_red = LinearSegmentedColormap('BLueRed', cdict)
MU = 4*math.pi*1e-7

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
    
class ConstantLoad(DessiaObject):
    """ 
    Sets a load on the selected elements by imposing a source value for the \
    current density vector J. Each element creates a source of the input value. 
    
    :param elements: The triangular elements.
    :type elements: List of volmdlr.TriangularElements objects
    :param value: Set the elements' current density vector J value.
    :type value: float
    """
    def __init__(self, elements: List[vmmesh.TriangularElement],
                 value: float):
        self.elements = elements
        self.value = value
        
        self.value_per_element = []
        total_area = sum([elem.area for elem in self.elements])
        for element in self.elements:
            self.value_per_element.append(value * element.area/total_area)
        DessiaObject.__init__(self, name='')
class SingleNodeLoad(DessiaObject):
    """ 
    Forces the value of the vector potential A at a node. To set a magnetic wall \
    the value of the vector potential has to be set to A.
    
    :param node: The node.
    :type node: volmdlr.Point2D object
    :param value: Set the node's vector potential A value.
    :type value: float
    """
    def __init__(self, node: vm.Point2D, value: float):
        self.node = node
        self.value = value
        
        DessiaObject.__init__(self, name='')
        
    
class MagnetLoad(DessiaObject):
    """
    Sets a load on the selected elements by imposing a source value for the \
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
             and (linear_element.points[0] not in self.non_contour_nodes \
             or linear_element.points[1] not in self.non_contour_nodes):
                contour_linear_elements.append(linear_element)
        return contour_linear_elements
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
class FiniteElementAnalysis(DessiaObject):
    """
    :param mesh: The meshing of the machine.
    :type mesh: Mesh object
    :param element_loads: The list of the loads applied to the triangluar elements.
    :type element_loads: List of ConstantLoad objects
    :param node_loads: The list of the loads applied to the nodes.
    :type node_loads: List of SingleNodeLoad objects
    :param magnet_loads: The list of the loads applied to the triangular elements.
    :type magnet_loads: List of MagnetLoad objects
    :param continuity_conditions: The list of continuity conditions applied to the nodes.
    :type continuity_conditions: List of ContinuityCondition objects
    """
    def __init__(self, mesh: vmmesh.Mesh,
                 element_loads: List[ConstantLoad],
                 node_loads: List[SingleNodeLoad],
                 magnet_loads: List[MagnetLoad],
                 continuity_conditions: List[ContinuityCondition]):
        self.mesh = mesh
        self.element_loads = element_loads  # current density J
        self.node_loads = node_loads 
        self.magnet_loads = magnet_loads
        self.continuity_conditions = continuity_conditions
        
        self.nb_loads = len(node_loads)
        
        DessiaObject.__init__(self, name='')
        
    def create_matrix(self):
        row_ind = []
        col_ind = []
        data = []
        
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
                
                row_ind.extend((indexes[0], indexes[0], indexes[0], indexes[1], indexes[1], indexes[1], indexes[2], indexes[2], indexes[2]))
                col_ind.extend((indexes[0], indexes[1], indexes[2], indexes[0], indexes[1], indexes[2], indexes[0], indexes[1], indexes[2]))
                data.extend((1/elements_group.mu_total * (b1**2 + c1**2) * element.area, 
                             1/elements_group.mu_total * (b1*b2 + c1*c2) * element.area,
                             1/elements_group.mu_total * (b1*b3 + c1*c3) * element.area,
                             1/elements_group.mu_total * (b1*b2 + c1*c2) * element.area,
                             1/elements_group.mu_total * (b2**2 + c2**2) * element.area,
                             1/elements_group.mu_total * (b2*b3 + c2*c3) * element.area,
                             1/elements_group.mu_total * (b1*b3 + c1*c3) * element.area,
                             1/elements_group.mu_total * (b2*b3 + c2*c3) * element.area,
                             1/elements_group.mu_total * (b3**2 + c3**2) * element.area))
                
        for i, load in enumerate(self.node_loads):
            index = self.mesh.node_to_index[load.node]
            
            row_ind.extend((len(self.mesh.nodes) + i, index))
            col_ind.extend((index, len(self.mesh.nodes) + i))
            data.extend((1, 1))
            
        for i, condition in enumerate(self.continuity_conditions):
            index1 = self.mesh.node_to_index[condition.node1]
            index2 = self.mesh.node_to_index[condition.node2]
            
            row_ind.extend((len(self.mesh.nodes)+len(self.node_loads)+i,
                            index1,
                            len(self.mesh.nodes)+len(self.node_loads)+i,
                            index2))
            col_ind.extend((index1,
                            len(self.mesh.nodes)+len(self.node_loads)+i,
                            index2,
                            len(self.mesh.nodes)+len(self.node_loads)+i))
            data.extend((1, 1, -condition.value, -condition.value))
            
        matrix = sparse.csr_matrix((data, (row_ind, col_ind)))
        return matrix
            
    def create_source_matrix(self):
        matrix = npy.zeros((len(self.mesh.nodes)+self.nb_loads+len(self.continuity_conditions), 1))
        for load in self.element_loads:
            for element in load.elements:
                indexes = [self.mesh.node_to_index[element.points[0]],
                           self.mesh.node_to_index[element.points[1]],
                           self.mesh.node_to_index[element.points[2]]]
                
                x1 = element.points[0][0]
                y1 = element.points[0][1]
                x2 = element.points[1][0]
                y2 = element.points[1][1]
                x3 = element.points[2][0]
                y3 = element.points[2][1]
                
                det_jacobien = abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
                
                element_form_functions = element.form_functions
                a1 = element_form_functions[0][0]
                b1 = element_form_functions[0][1]
                c1 = element_form_functions[0][2]
                a2 = element_form_functions[1][0]
                b2 = element_form_functions[1][1]
                c2 = element_form_functions[1][2]
                a3 = element_form_functions[2][0]
                b3 = element_form_functions[2][1]
                c3 = element_form_functions[2][2]
                
                double_integral_N1_dS = det_jacobien*(a1 + 0.5*b1*x2 + 0.5*c1*y2 + 0.5*b1*x3 + 0.5*c1*y3)
                double_integral_N2_dS = det_jacobien*(a2 + 0.5*b2*x2 + 0.5*c2*y2 + 0.5*b2*x3 + 0.5*c2*y3)
                double_integral_N3_dS = det_jacobien*(a3 + 0.5*b3*x2 + 0.5*c3*y2 + 0.5*b3*x3 + 0.5*c3*y3)
                
                matrix[indexes[0]][0] += load.value * double_integral_N1_dS 
                matrix[indexes[1]][0] += load.value * double_integral_N2_dS 
                matrix[indexes[2]][0] += load.value * double_integral_N3_dS
                
        for i, load in enumerate(self.node_loads):
            matrix[len(self.mesh.nodes) + i][0] += load.value
            
        for magnet_load in self.magnet_loads:
            for linear_element in magnet_load.contour_linear_elements():
                indexes = [self.mesh.node_to_index[linear_element.points[0]],
                           self.mesh.node_to_index[linear_element.points[1]]]
                length = linear_element.length()
                dl = vm.Vector2D([-linear_element.interior_normal[1],
                                  linear_element.interior_normal[0]])
                matrix[indexes[0]][0] += magnet_load.magnetization_vector.Dot(dl) * length/2
                matrix[indexes[1]][0] += magnet_load.magnetization_vector.Dot(dl) * length/2
                
        return matrix
    
    def solve(self):
        """
        Solve the matix equation : F = K.X, where X is the unknown vector. \
        Each value of the vector represent the amplitude of the vector potential \
        of a node. 
        Returns a Result object.
        """
        # TEST RESULT : 
        # Between :
        # - npy.linalg.solve
        # - scipy.linalg.solve
        # - scipy.linalg.solve(assume_a='sym')
        # - scipy.sparse.linalg.spsolved
        # - scipy.pinvh than dot
        
        # scipy.sparse.linalg.spsolved is the fastest !
        # print('avant')
        K_sparse = self.create_matrix()
        F = self.create_source_matrix()
        try:
            X = sparse.linalg.spsolve(K_sparse, F,
                                      permc_spec='NATURAL',
                                      use_umfpack=True)
        except sparse.linalg.MatrixRankWarning:
            print('MatricRankWarning')
            raise NotImplementedError
        X = list(X)
        # print('apres')
        return Result(self.mesh, X)
        
    def plot_elements_loads(self, ax=None):
        """ 
        Plots the mesh. The triangular elements are filled red if they are a \
        source of magnetic potential thanks to their current density.
        """
        if ax is None:
            fig, ax = plt.subplots()
        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:
                element.plot(ax=ax, color='w', fill=True)
        
        for load in self.element_loads:
            for element in load.elements:
                element.plot(ax=ax, color='r', fill=True)
        return ax
    
    def plot_elements_permeability(self, ax=None):
        """ 
        Plots the mesh with colored triangular elements depending on the \ 
        permeability of the material the element meshes. 
        """
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()
        
        color_map = ((0, 0, 1), (1, 0, 0))
        permeabilities = []
        for elements_group in self.mesh.elements_groups:
            permeabilities.append(elements_group.mu_total)
        mu_max = max(permeabilities)
        mu_min = min(permeabilities)
        colors = []
        for elements_group in self.mesh.elements_groups:
            x = (elements_group.mu_total - mu_min) / (mu_max - mu_min)
            color = (color_map[0][0]-(color_map[0][0]-color_map[1][0])*x, 
                     color_map[0][1]-(color_map[0][1]-color_map[1][1])*x,
                     color_map[0][2]-(color_map[0][2]-color_map[1][2])*x)
            colors.append(color)
        
        norm = mpl.colors.Normalize(vmin=mu_min, vmax=mu_max)
        sm = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)  #plt.get_cmap('jet_r')
        sm.set_array([])
        cbar = fig.colorbar(sm, ticks=npy.linspace(mu_min, mu_max, 10))
        cbar.set_label('permeability')
        for i, elements_group in enumerate(self.mesh.elements_groups):
            for element in elements_group.elements:
                element.plot(ax=ax, color=colors[i], fill=True)
                
        return ax
    
    def plot_magnet_loads(self, ax=None):
        """ 
        Plots the mesh. The triangular elements are filled red if they are a \
        source of magnetic potential thanks to their magnetization. The contour \
        of the magnetizating mesh is drawn in blue. 
        """
        if ax is None:
            fig, ax = plt.subplots()
        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:
                element.plot(ax=ax, color='w', fill=True)
            
        for magnet_load in self.magnet_loads:
            for element in magnet_load.elements:
                element.plot(ax=ax, color='r', fill=True)
                
            contour_linear_elements = magnet_load.contour_linear_elements()
            for linear_element in contour_linear_elements:
                linear_element.plot(ax=ax, color='b')
        
        print(self.magnet_loads)
        
        return ax
    
    def plot_continuity_condition(self, ax=None):
        if ax is None:
            ax = self.mesh.plot()
            
        for i, continuity_condition in enumerate(self.continuity_conditions):
            continuity_condition.node1.MPLPlot(ax=ax, color='C{}'.format(i % 10))
            continuity_condition.node2.MPLPlot(ax=ax, color='C{}'.format(i % 10))
        return ax
                
    
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
    def plot_brbtetha(self, ax=None, air_gap_elements_group_name='Gap ring'):
        if ax is None:
            fig, ax = plt.subplots()
            # ax = self.mesh.plot()
        else:
            fig = plt.gcf()
            # self.mesh.plot(ax=ax)
        element_to_magnetic_field = self.magnetic_field_per_element()
        
        for elements_group in self.mesh.elements_groups:
            if elements_group.name == air_gap_elements_group_name:
                gap_elements_group = elements_group
                break
        
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
            
        color_map = ((0,0,1), (1,0,0))
        Bs = [B.norm() for B in list(element_to_magnetic_field.values())]
        
        if Bmax is None and Bmin is None:
            B_max, B_min = max(Bs), min(Bs)
        elif Bmax is not None and Bmin is None:
            B_max, B_min = Bmax, min(Bs)
        elif Bmax is None and Bmin is not None:
            B_max, B_min = max(Bs), Bmin
        else:
            B_max, B_min = Bmax, Bmin
            
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
        
        color_map = ((0,0,1), (1,0,0))
        Bs = [B.norm() for B in list(element_to_magnetic_field.values())]
        
        if Bmax is None:
            B_max, B_min = max(Bs), min(Bs)
        else:
            B_max, B_min = Bmax, min(Bs)
        
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
        
        from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner
        import matplotlib.cm as cm
            
        element_to_magnetic_field = self.magnetic_field_per_element()
        
        x = []
        y = []
        Z = []
        for group in self.mesh.elements_groups:
            for element in group.elements:
                x_center, y_center = element.center
                x.append(x_center)
                y.append(y_center)
                Z.append(element_to_magnetic_field[element].Norm())
        
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
        x_list = []
        y_list = []
        for node in self.mesh.nodes:
            x_list.append(node[0])
            y_list.append(node[1])
        triangles = []
        for group in self.mesh.elements_groups:
            for element in group.elements:
                triangles.append([self.mesh.node_to_index[element.points[0]], 
                                  self.mesh.node_to_index[element.points[1]],
                                  self.mesh.node_to_index[element.points[2]]])
            
        x = npy.asarray(x_list)
        y = npy.asarray(y_list)
        z = npy.asarray([p[0] for p in self.result_vector[:len(self.mesh.nodes)]])
        z_min, z_max = min(z), max(z)
        triang = mtri.Triangulation(x, y, triangles)
        
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