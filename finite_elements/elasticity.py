"""
Created on Mon Aug 31 15:45:15 2020

@author: gasmi

"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as npy
import matplotlib.tri as mtri
import volmdlr as vm
import volmdlr.mesh as vmmesh
import finite_elements.core as corefe
import math
from scipy import sparse
from scipy import linalg
import time  
from dessia_common import DessiaObject
from typing import TypeVar, List, Tuple,Dict
from itertools import product

cdict = {'red':  [(0.0, 0.0, 0.0),
                   (1.0, 1.0, 1.0)],
         'green': [(0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)],
         'blue':  [(0.0, 1.0, 1.0),
                   (1.0, 0.0, 0.0)]}
blue_red = LinearSegmentedColormap('BLueRed', cdict)

    
class ConstantLoad(DessiaObject):
    """ 
    Sets a load on the selected elements by imposing a source value for the \
    current load vector F . Each element creates a source of the input value. 
    
    :param elements: The triangular elements.
    :type elements: List of volmdlr.TriangularElements objects
    :param values: Set the elements' current load vector F value along both axis.
    :type values: List of float
    """
    def __init__(self, elements:List[vmmesh.TriangularElement], values:List[float]):
        self.elements = elements
        self.values = values
        
        self.value_per_element = []
        total_area = sum([elem.area for elem in self.elements])
        for element in self.elements:
            self.value_per_element.append(values[0] * element.area/total_area)
            self.value_per_element.append(values[1] * element.area/total_area)

        DessiaObject.__init__(self, name='')
        
        
        
class BoundaryLoad(DessiaObject):
        """ 
    Sets a load by unit of length on the selected bound (a linear element) of the domain.

    :param point1: The first point of the linear element.
    :type point1: vm.Point2D object
    :param point2: The first point of the linear element.
    :type point2: vm.Point2D object
    :param interior_normal: The interior normal of the linear element
    :type interior_normal: vm.Vector2D object
    :param load_vector: Value of the force in both directions
    :type load_vector: vm.Vector2D object
    """
        def __init__(self,point1:vm.Point2D,point2:vm.Point2D,interior_normal:vm.Vector2D,load_vector:vm.Vector2D):
            self.point1=point1
            self.point2=point2
            self.interior_normal=interior_normal
            self.linear_element=vmmesh.LinearElement([self.point1,self.point2],self.interior_normal)
            self.load_vector=load_vector
      
class SingleNodeLoad(DessiaObject):
    """
    Forces the value of the vector potential A at a node. To set a magnetic wall \
    the value of the vector potential has to be set to A.
    
    :param node: The node.
    :type node: volmdlr.Point2D object
    :param value: Set the node's vector potential A value.
    :type value: float
    """
    def __init__(self, node:vm.Point2D, values:List[float]):
        self.node = node
        self.values = values
        
        DessiaObject.__init__(self, name='')
        
class NodeDisplacement(DessiaObject):
    """ 
    Forces the value of the displacement u along the x axis and v along the/
    y axis at a node. 

    :param node: The node.
    :type node: volmdlr.Point2D object
    :param u: Set the node's displacement along the x axis.
    :type value: float
    :param v: Set the node's displacement along the y axis.
    :type value: float
    """
    def __init__(self, node:vm.Point2D, u:float,v:float):
        self.node = node
        self.u = u
        self.v=v
        
        DessiaObject.__init__(self, name='')
        


class FiniteElementAnalysis(DessiaObject):
    """
    :param mesh: The meshing of the machine.
    :type mesh: Mesh object
    :param materials: Dictionnary of elements groups and their respective young module and poisson coefficient 
    :type materials:  Dictionnary of elements groups and a list of their young modules and poisson coefficients
    :param element_loads: The list of the loads applied to the triangluar elements.
    :type element_loads: List of ConstantLoad objects
    :param node_loads: The list of the loads applied to the nodes.
    :type node_loads: List of SingleNodeLoad objects
    :param continuity_conditions: The list of continuity conditions applied to the nodes.
    :type continuity_conditions: List of ContinuityCondition objects
    """
    def __init__(self,mesh:vmmesh.Mesh,materials:corefe.Materials,element_loads:List[ConstantLoad],boundary_loads:[BoundaryLoad],node_loads:List[SingleNodeLoad],node_displacements:List[NodeDisplacement]):
        self.mesh = mesh
        self.materials=materials
        self.node_loads=node_loads
        self.element_loads = element_loads 
        self.boundary_loads=boundary_loads
        self.node_displacements= node_displacements
        self.nb_loads = len(node_loads)
  
        
        DessiaObject.__init__(self, name='')
        
    def create_stiffness_matrix(self):
        """
        Creats the elementary stiffness matrix Ke of each element and assembles/
        them into the matrix K. The rows and columns corresponding to the imposed/
        node displacements are also added (lagrange multiplicators method).
       
        """
        row_ind=[]
        col_ind=[]
        data = []       
               
        for elements_group in self.mesh.elements_groups:
            e_young=self.materials.materials_properties[elements_group][0]
            v_poisson=self.materials.materials_properties[elements_group][1]
            alpha=e_young/(1-v_poisson**2)
            beta=(1-v_poisson)/2
            gamma=(1+v_poisson)/2
            for element in elements_group.elements:
                indexes=[]
                element_form_functions = element.form_functions
                for point in element.points:
            
                  indexes.append(self.mesh.set_node_displacement_index()[point][0])
                  indexes.append(self.mesh.set_node_displacement_index()[point][1])
                   # x=2*self.mesh.node_to_index[point]
                   # y=2*self.mesh.node_to_index[point]+1
                   # indexes.append(x)
                   # indexes.append(y)
               
                
                b1 = element_form_functions[0][1] 
                c1 = element_form_functions[0][2]
                b2 = element_form_functions[1][1]
                c2 = element_form_functions[1][2]
                b3 = element_form_functions[2][1]
                c3 = element_form_functions[2][2]
               
                
                
                for k in range(len(indexes)):
                    
                    col_ind.extend((indexes[k])for k in range(len(indexes)))
                    
                for k,u in enumerate(list(product(indexes,repeat=2))):
                    
                    row_ind.append(u[0])
                
                
               
                data.extend((element.area*(b1**2+beta*c1**2)*alpha,
                             element.area*gamma*b1*c1*alpha,
                             element.area*(b1*b2+beta*c2*c1)*alpha,
                             element.area*(c2*b1*v_poisson+beta*b2*c1)*alpha,
                             element.area*(b1*b3+beta*c3*c1)*alpha,
                             element.area*(v_poisson*c3*b1+beta*b3*c1)*alpha,
                             
                             element.area*gamma*b1*c1*alpha,
                             element.area*(c1**2+beta*b1**2)*alpha,
                             element.area*(v_poisson*c1*b2+beta*c2*b1)*alpha,
                             element.area*(c2*c1+beta*b2*b1)*alpha,
                             element.area*(v_poisson*c1*b3+beta*c3*b1)*alpha,
                             element.area*(c1*c3+beta*b3*b1)*alpha,
                             
                             
                             element.area*(b1*b2+beta*c1*c2)*alpha,
                             element.area*(v_poisson*c1*b2+beta*b1*c2)*alpha,
                             element.area*(b2**2+beta*c2**2)*alpha,
                             element.area*b2*c2*gamma*alpha,
                             element.area*(b3*b2+beta*c3*c2)*alpha,  
                             element.area*(v_poisson*c3*b2+beta*c2*b3)*alpha,   
                             
                             element.area*(c2*b1*v_poisson+beta*b2*c1)*alpha,
                             element.area*(c2*c1+beta*b2*b1)*alpha,
                             element.area*b2*c2*gamma*alpha,
                             element.area*(c2**2+beta*b2**2)*alpha,
                             element.area*(v_poisson*b3*c2+beta*c3*b2)*alpha,  
                             element.area*(c3*c2+beta*b3*b2)*alpha,   
                             
                             element.area*(b1*b3+beta*c1*c3)*alpha,
                             element.area*(v_poisson*c1*b3+beta*b1*c3)*alpha,
                             element.area*(b3*b2+beta*c3*c2)*alpha,
                             element.area*(v_poisson*b3*c2+beta*c3*b2)*alpha, 
                             element.area*(b3**2+beta*c3**2)*alpha,
                             element.area*b3*c3*gamma*alpha,
                             
                             element.area*(v_poisson*b1*c3+beta*b3*c1)*alpha,
                             element.area*(beta*b1*b3+c1*c3)*alpha,
                             element.area*(v_poisson*c3*b2+beta*b3*c2)*alpha,
                             element.area*(c2*c3+beta*b2*b3)*alpha, 
                             element.area*b3*c3*gamma**alpha,
                             element.area*(c3**2+beta*b3**2)*alpha))
                             
                             
               
                
            
               

        
        for  i,displacement in enumerate(self.node_displacements):
            
                  k = self.node_displacements.index(displacement)
                  if displacement.u == None and displacement.v == None :
                      break 
                  if displacement.u!=None:
                      index = self.mesh.node_to_index[displacement.node]
                    
                
            
                      row_ind.extend((2*len(self.mesh.nodes) + i+k, 2*index,2*len(self.mesh.nodes)+i+k,2*index))
                      col_ind.extend((2*index, 2*len(self.mesh.nodes) + i+k,2*index,2*len(self.mesh.nodes)+i+k))
            
                      data.extend((1,1,0,0))
                
                  if displacement.v!=None:
                    index = self.mesh.node_to_index[displacement.node]
                
                
            
                    row_ind.extend((2*len(self.mesh.nodes) + i+1+k, 2*index+1,2*len(self.mesh.nodes)+i+1+k,2*index+1))
                    col_ind.extend((2*index+1, 2*len(self.mesh.nodes) + i+1+k,2*index+1,2*len(self.mesh.nodes)+i+1+k))
            
                    data.extend((1,1,0,0))
 
            
        matrix = sparse.csr_matrix((data, (row_ind, col_ind)))  
      
        

        return matrix
         
    
    def create_force_matrix(self):
        
        matrix = npy.zeros((2*len(self.mesh.nodes)+2*len(self.node_displacements), 1))
        
        for load in self.element_loads:
                for element in load.elements:
                    indexes=[]
                    for point in element.points:
                        indexes.append(self.mesh.set_node_displacement_index()[point][0])
                        indexes.append(self.mesh.set_node_displacement_index()[point][1])
                  
                
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
                
                I=[double_integral_N1_dS,double_integral_N2_dS,double_integral_N3_dS]
                k=0
                i=0
                while k <len(indexes): 
             
                        matrix[indexes[k]][0] += load.values[0] * I[i]
                        matrix[indexes[k]][0] += load.values[1] * I[i]
                        i=i+1
                        k=k+2
                       
                          
                
        for load in self.node_loads:
            indexes=[]
      
            
                  
            x=2*self.mesh.node_to_index[load.node]
            y=2*self.mesh.node_to_index[load.node]+1
                   
            indexes.append(x)
            indexes.append(y)  
            
            
            matrix[x][0] += load.values[0]
            matrix[y][0]+=load.values[1]
            
            
            
        for  i,displacement in enumerate(self.node_displacements):
             k=self.node_displacements.index(displacement)
            
             matrix[ 2*len(self.mesh.nodes) + i+k][0]=displacement.u
             matrix[ 2*len(self.mesh.nodes) + i+1+k][0]=displacement.v
             
             
        for boundary_load in self.boundary_loads:
           
                indexes =[]
                indexes.append(self.mesh.set_node_displacement_index()[point][0])
                indexes.append(self.mesh.set_node_displacement_index()[point][1])
               
                length = boundary_load.linear_element.length()
                dl = vm.Vector2D([boundary_load.interior_normal[0],boundary_load.interior_normal[1]])
                matrix[indexes[0]][0] += boundary_load.load_vector.Dot(dl) * length/2
                matrix[indexes[1]][0] += boundary_load.load_vector.Dot(dl) * length/2
                
        return matrix   
    
    def solve(self):
        """
        Solve the matix equation : F = K.X, where X is the unknown vector. \
        Each pair of value of the vector represents the displacement of a node\
        along the x and y axis.
        
        Returns a Result object.
        """
        
        K_sparse = self.create_stiffness_matrix()
        F = self.create_force_matrix()
        
        try:
            X = sparse.linalg.spsolve(K_sparse, F, permc_spec='NATURAL',use_umfpack=True)
           
        except sparse.linalg.MatrixRankWarning:
            print('MatricRankWarning')
            raise NotImplementedError
        X = list(X)
        
        return  FiniteElementAnalysisResult(self.mesh,self.materials, X)
        
    def plot_elements_loads(self, ax=None):
        """ 
        Plots the mesh. The triangular elements are filled red if they are a \
       loaded
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
  
    
    
        
class  FiniteElementAnalysisResult(DessiaObject):
    """
    :param mesh: The mesh on which the result has been solved.
    :type mesh: a Mesh object
    :param result_vector: The solution vector for the displacement U.
    :type result_vector: a list of float
    """
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True
    def __init__(self, mesh:vmmesh.Mesh,materials:corefe.Materials,result_vector:List[float]):
        self.mesh = mesh
        self.materials=materials
        self.result_vector = result_vector
        
        DessiaObject.__init__(self, name='')
        
    def deformation_field_per_element(self):
        element_to_deformation_field = {}
        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:
                element_form_functions = element.form_functions
    
                
                b1 = element_form_functions[0][1]
                c1 = element_form_functions[0][2]
                b2 = element_form_functions[1][1]
                c2 = element_form_functions[1][2]
                b3 = element_form_functions[2][1]
                c3 = element_form_functions[2][2]            
               
                gradient_element=[[0,0.5*b1,0.5*c1],
                                  [0,0.5*b2,0.5*c2],
                                  [0,0.5*b3,0.5*c3]]
                                                                   
                                                   
                element_to_deformation_field[element]=(gradient_element+
                                                       npy.transpose(gradient_element))
                
                                                                                                                                                                                                                                                                                                                                                           
        return element_to_deformation_field
  
    def von_mises_stress_per_element(self):
        "Von mises stress to represent the stress tensor"                                        
        element_to_von_mises_stress ={}
        for elements_group in self.mesh.elements_groups:
            e_young=self.materials.materials_properties[elements_group][0]
            v_poisson=self.materials.materials_properties[elements_group][1]
            mu_lame=e_young/(2*(1+v_poisson))
            for element in elements_group.elements:
                
              
                element_to_deformation_field=self.deformation_field_per_element()
                
                sigma=[2*s*mu_lame for s in element_to_deformation_field[element]]
                w=linalg.eigh(sigma,eigvals_only=True) 
                vm_stress=(1/math.sqrt(2))*math.sqrt((w[0]-w[1])**2+(w[1]-w[2])**2+
                                                     (w[0]-w[2])**2) 
                element_to_von_mises_stress[element]=vm_stress
                
        return  element_to_von_mises_stress
               
    
    def plot_stress_field(self, ax=None, stress_max=None): 
        """
       Plots the mesh with colored triangular elements representing the \
    value of the Von mises stress inside the structure. 
      """
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()

        element_to_von_mises_stress = self.von_mises_stress_per_element()
        
        color_map = ((0,0,1), (1,0,0))
        vm_stresses = [vm_stress for vm_stress in list(element_to_von_mises_stress.values())]
        
        if stress_max is None:
            stress_max, stress_min = max(vm_stresses), min(vm_stresses)
        else:
            stress_max, stress_min = stress_max, min(vm_stresses)
        
        stress_to_color = {}
        for stress in vm_stresses:
            if stress > stress_max:
                x = 1
            else:
                x = (stress - stress_min) / (stress_max - stress_min)
            color = (color_map[0][0]-(color_map[0][0]-color_map[1][0])*x, 
                     color_map[0][1]-(color_map[0][1]-color_map[1][1])*x,
                     color_map[0][2]-(color_map[0][2]-color_map[1][2])*x)
            stress_to_color[stress] = color
        
        for group in self.mesh.elements_groups:
            for element in group.elements:
                element.plot(ax=ax, color=stress_to_color[element_to_von_mises_stress[element]],fill=True)
        
        norm = mpl.colors.Normalize(vmin=stress_min, vmax=stress_max)
        sm = plt.cm.ScalarMappable(cmap=blue_red, norm=norm) 
        sm.set_array([])
        cbar = fig.colorbar(sm, ticks=npy.linspace(stress_min, stress_max, 10))
        cbar.set_label('Von mises Stress in Pa')
        
        return ax
    
    def node_to_displacement(self):
        """
        Creats a dictionnary containing : for keys the nodes of the mesh and for/
        values their respective displacement along the x and y axes.
        """
        U_coordonate={}
        nodes = self.mesh.nodes
       
        for node in nodes:
             indexes=[]
             indexes.append(self.mesh.set_node_displacement_index()[node][0])
             indexes.append(self.mesh.set_node_displacement_index()[node][1])
                            
             U_coordonate[node]=[self.result_vector[indexes[0]],self.result_vector[indexes[1]]]
                        
        return U_coordonate

                        

       