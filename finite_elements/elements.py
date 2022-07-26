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
from dessia_common import DessiaObject
# from typing import List #Tuple, TypeVar
import numpy as npy
import finite_elements.core


class MagneticElement2D(vmmesh.TriangularElement2D):
    """
    """
    # _standalone_in_db = False
    # _non_serializable_attributes = []
    # _non_eq_attributes = ['name']
    # _non_hash_attributes = ['name']
    # _generic_eq = True
    def __init__(self, triangular_element: vmmesh.TriangularElement2D,
                 mu_total: float, name : str = ''):
        self.triangular_element = triangular_element
        vmmesh.TriangularElement2D.__init__(self, points=triangular_element.points)
        self.mu_total = mu_total

        # DessiaObject.__init__(self, name=name)

    @property
    def dimension(self):
        return 1

    def elementary_matrix(self):
        """
        Create the elementary matrix of the MagneticElement2D

        :return: data
        """

        element_form_functions = self.triangular_element.form_functions
        b1 = element_form_functions[0][1]
        c1 = element_form_functions[0][2]
        b2 = element_form_functions[1][1]
        c2 = element_form_functions[1][2]
        b3 = element_form_functions[2][1]
        c3 = element_form_functions[2][2]

        data = (1/self.mu_total * (b1**2 + c1**2) * self.triangular_element.area,
                1/self.mu_total * (b1*b2 + c1*c2) * self.triangular_element.area,
                1/self.mu_total * (b1*b3 + c1*c3) * self.triangular_element.area,
                1/self.mu_total * (b1*b2 + c1*c2) * self.triangular_element.area,
                1/self.mu_total * (b2**2 + c2**2) * self.triangular_element.area,
                1/self.mu_total * (b2*b3 + c2*c3) * self.triangular_element.area,
                1/self.mu_total * (b1*b3 + c1*c3) * self.triangular_element.area,
                1/self.mu_total * (b2*b3 + c2*c3) * self.triangular_element.area,
                1/self.mu_total * (b3**2 + c3**2) * self.triangular_element.area)

        return data

    def elementary_source_matrix(self, indexes):
        """
        Create the elementary source matrix of the MagneticElement2D

        :return: (double_integral_N1_dS, double_integral_N2_dS, double_integral_N3_dS)
        """

        x1 = self.triangular_element.points[0][0]
        y1 = self.triangular_element.points[0][1]
        x2 = self.triangular_element.points[1][0]
        y2 = self.triangular_element.points[1][1]
        x3 = self.triangular_element.points[2][0]
        y3 = self.triangular_element.points[2][1]

        det_jacobien = abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))

        element_form_functions = self.triangular_element.form_functions
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

        return (double_integral_N1_dS, double_integral_N2_dS, double_integral_N3_dS)


# class MagneticElementsGroup(vmmesh.ElementsGroup):
#     # _standalone_in_db = False
#     # _non_serializable_attributes = []
#     # _non_eq_attributes = ['name']
#     # _non_hash_attributes = ['name']
#     # _generic_eq = True
#     def __init__(self, magnetic_elements: List[MagneticElement2D],
#                  mu_total: float, name: str):
#         self.magnetic_elements = magnetic_elements
#         self.triangular_elements = self._triangular_elements()
#         vmmesh.ElementsGroup.__init__(self, elements=self.triangular_elements, name=name)
#         self.mu_total = mu_total

#         # DessiaObject.__init__(self, name=name)

#     def _triangular_elements(self):
#         return [element.triangular_element for element in self.magnetic_elements]


class ElasticityElement(DessiaObject):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True
    # def __init__(self, mesh_element,
    #              elasticity_modulus : float,
    #              poisson_ratio : float,
    #              name : str = ''):
    #     self.mesh_element = mesh_element
    #     self.elasticity_modulus = elasticity_modulus
    #     self.poisson_ratio = poisson_ratio

    def __init__(self, mesh_element,
                 elasticity_modulus, poisson_ratio,
                 mass_density,
                 displacements = None,
                 stress = None,
                 strain = None,
                 name : str = ''):
        self.mesh_element = mesh_element
        self.elasticity_modulus = elasticity_modulus
        self.poisson_ratio = poisson_ratio
        self.mass_density = mass_density
        self.points = self._points()
        self.b_matrix = self._b_matrix()
        self.d_matrix_plane_strain = self._d_matrix_plane_strain()
        self.d_matrix_plane_stress = self._d_matrix_plane_stress()
        self.displacements = displacements
        self.stress = stress
        self.strain = strain

        DessiaObject.__init__(self, name=name)

    def d_matrix(self, plane_strain: bool, plane_stress:bool):

        if (plane_strain and plane_stress):
            raise ValueError('just one of plane_strain or plane_stress can be True')
        elif (not plane_strain and not plane_stress):
            raise ValueError('one of plane_strain or plane_stress must be True')
        elif plane_strain:
            return self.d_matrix_plane_strain
        elif plane_stress:
            return self.d_matrix_plane_stress

    def _points(self):
        obj=getattr(finite_elements.core, f'Node{self.__class__.__name__[-2::]}')
        return [getattr(obj, 'from_point')(point) for point in self.mesh_element.points]

class ElasticityTriangularElement2D(ElasticityElement, vmmesh.TriangularElement2D):
    # _standalone_in_db = False
    # _non_serializable_attributes = []
    # _non_eq_attributes = ['name']
    # _non_hash_attributes = ['name']
    # _generic_eq = True
    # def __init__(self, mesh_element: vmmesh.TriangularElement2D,
    #              elasticity_modulus, poisson_ratio,
    #              thickness: float = 1.0,
    #              displacements = None,
    #              stress = None,
    #              strain = None,
    #              name : str = ''):
    #     self.mesh_element = mesh_element
    #     self.elasticity_modulus = elasticity_modulus
    #     self.poisson_ratio = poisson_ratio
    #     self.thickness = thickness
    #     self.points = self.mesh_element.points
    #     self.b_matrix = self._b_matrix()
    #     self.d_matrix = self._d_matrix()
    #     self.displacements = displacements
    #     self.stress = stress
    #     self.strain = strain

    def __init__(self, mesh_element: vmmesh.TetrahedralElement,
                 elasticity_modulus, poisson_ratio, mass_density,
                 thickness: float = 1.0,
                 displacements = None,
                 stress = None,
                 strain = None,
                 name : str = ''):
        self.thickness = thickness

        ElasticityElement.__init__(self, mesh_element,
                                       elasticity_modulus, poisson_ratio, mass_density)
        vmmesh.TriangularElement2D.__init__(self, points=mesh_element.points)

        # DessiaObject.__init__(self, name=name)

    def _b_matrix(self):

        y = [(self.points[i].y-self.points[j].y) for (i,j) in [(1,2), (2,0), (0,1)]]
        x = [(self.points[i].x-self.points[j].x) for (i,j) in [(2,1), (0,2), (1,0)]]

        det_jacobian = (self.points[0].x-self.points[2].x)*(self.points[1].y-self.points[2].y) \
            - (self.points[0].y-self.points[2].y)*(self.points[1].x-self.points[2].x)

        data = [y[0], 0, y[1], 0, y[2], 0,
                0, x[0], 0, x[1], 0, x[2],
                x[0], y[0], x[1], y[1], x[2], y[2]]

        b_matrix = (1/det_jacobian) * npy.array(data).reshape(3,6)

        return b_matrix

    def _d_matrix_plane_strain(self):

        a = (self.elasticity_modulus*self.poisson_ratio) \
            /((1+self.poisson_ratio)*(1-2*self.poisson_ratio))
        b = self.elasticity_modulus/(2*(1+self.poisson_ratio))
        c = a+2*b
        data = [c, a, 0,
                a, c, 0,
                0, 0, b]

        return npy.array(data).reshape(3,3)

    def _d_matrix_plane_stress(self):

        a = (self.elasticity_modulus*self.poisson_ratio) \
            /(1-self.poisson_ratio**2)
        b = self.elasticity_modulus/(2*(1+self.poisson_ratio))
        c = a+2*b

        data = [c, a, 0,
                a, c, 0,
                0, 0, b]

        return npy.array(data).reshape(3,3)

    @property
    def dimension(self):
        return 2

    def elementary_matrix(self, plane_strain: bool, plane_stress:bool):

        b_matrix = self.b_matrix
        d_matrix = self.d_matrix(plane_strain, plane_stress)

        # if plane_strain:
        #     d_matrix = self.self.d_matrix_plane_strain
        # elif plane_stress:
        #     d_matrix = self.self.d_matrix_plane_stress

        # y = [(self.points[i].y-self.points[j].y) for (i,j) in [(1,2), (2,0), (0,1)]]
        # x = [(self.points[i].x-self.points[j].x) for (i,j) in [(2,1), (0,2), (1,0)]]

        # det_jacobian = (self.points[0].x-self.points[2].x)*(self.points[1].y-self.points[2].y) \
        #     - (self.points[0].y-self.points[2].y)*(self.points[1].x-self.points[2].x)

        # data = [y[0], 0, y[1], 0, y[2], 0,
        #         0, x[0], 0, x[1], 0, x[2],
        #         x[0], y[0], x[1], y[1], x[2], y[2]]

        # b_matrix = (1/det_jacobian) * npy.array(data).reshape(3,6)

        # elasticity_modulus = self.elasticity_modulus
        # poisson_ratio = self.poisson_ratio

        # data = [1, poisson_ratio, 0,
        #         poisson_ratio, 1, 0,
        #         0, 0, (1-poisson_ratio)/2]

        # d_matrix = (elasticity_modulus/(1 - (poisson_ratio)**2)) * npy.array(data).reshape(3,3)

        stiffness_matrix = self.thickness * self.area * (
            npy.matmul(npy.matmul(b_matrix.transpose(), d_matrix), b_matrix))

        return stiffness_matrix.flatten()

    def elementary_mass_matrix(self):
        data = [2, 0, 1, 0, 1, 0,
                0, 2, 0, 1, 0, 1,
                1, 0, 2, 0, 1, 0,
                0, 1, 0, 2, 0, 1,
                1, 0, 1, 0, 2, 0,
                0, 1, 0, 1, 0, 2]

        x1 = self.mesh_element.points[0][0]
        y1 = self.mesh_element.points[0][1]
        x2 = self.mesh_element.points[1][0]
        y2 = self.mesh_element.points[1][1]
        x3 = self.mesh_element.points[2][0]
        y3 = self.mesh_element.points[2][1]

        det_jacobien = (abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)))

        mass_matrix = det_jacobien * ((self.mass_density * self.area \
                        * self.thickness)/12) * npy.array(data).reshape(6, 6)

        # mass_matrix = 0.5* det_jacobien * ((self.mass_density * self.area \
        #                 * self.thickness)/12) * npy.array(data).reshape(6, 6)

        return mass_matrix.flatten()

    # def elementary_mass_matrix(self):
    #     data = [2, 0, 1, 0, 1, 0,
    #             0, 2, 0, 1, 0, 1,
    #             1, 0, 2, 0, 1, 0,
    #             0, 1, 0, 2, 0, 1,
    #             1, 0, 1, 0, 2, 0,
    #             0, 1, 0, 1, 0, 2]

    #     x1 = self.mesh_element.points[0][0]
    #     y1 = self.mesh_element.points[0][1]
    #     x2 = self.mesh_element.points[1][0]
    #     y2 = self.mesh_element.points[1][1]
    #     x3 = self.mesh_element.points[2][0]
    #     y3 = self.mesh_element.points[2][1]

    #     d = [(x2*y3 - x3*y2), (x3*y1 - x1*y3), (x1*y2 - x2*y1), 0, 0, 0,
    #          (y2-y3), (y3-y1), (y1-y2), 0, 0, 0,
    #          (x3-x2), (x1-x3), (x2-x1), 0, 0, 0,
    #          0, 0, 0, (x2*y3 - x3*y2), (x3*y1 - x1*y3), (x1*y2 - x2*y1),
    #          0, 0, 0, (y2-y3), (y3-y1), (y1-y2),
    #          0, 0, 0, (x3-x2), (x1-x3), (x2-x1)]

    #     A = 1/(2*self.area) * npy.array(d).reshape(6, 6)
    #     from numpy.linalg import inv
    #     A = inv(A)

    #     M = npy.array(data).reshape(6, 6)

    #     mass_matrix = ((self.mass_density * self.area \
    #                     * self.thickness)/12) * npy.matmul(A, M)

    #     mass_matrix = ((self.mass_density * self.area \
    #                     * self.thickness)/12) * (npy.matmul(npy.matmul(A.transpose(), M), A))

    #     return mass_matrix.flatten()


    # def elementary_mass_matrix(self):

    #     x1 = self.mesh_element.points[0][0]
    #     y1 = self.mesh_element.points[0][1]
    #     x2 = self.mesh_element.points[1][0]
    #     y2 = self.mesh_element.points[1][1]
    #     x3 = self.mesh_element.points[2][0]
    #     y3 = self.mesh_element.points[2][1]

    #     det_jacobien = abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))

    #     element_form_functions = self.mesh_element.form_functions
    #     a1 = element_form_functions[0][0]
    #     b1 = element_form_functions[0][1]
    #     c1 = element_form_functions[0][2]
    #     a2 = element_form_functions[1][0]
    #     b2 = element_form_functions[1][1]
    #     c2 = element_form_functions[1][2]
    #     a3 = element_form_functions[2][0]
    #     b3 = element_form_functions[2][1]
    #     c3 = element_form_functions[2][2]

    #     N1 = det_jacobien*(a1 + 0.5*b1*x2 + 0.5*c1*y2 + 0.5*b1*x3 + 0.5*c1*y3)
    #     N2 = det_jacobien*(a2 + 0.5*b2*x2 + 0.5*c2*y2 + 0.5*b2*x3 + 0.5*c2*y3)
    #     N3 = det_jacobien*(a3 + 0.5*b3*x2 + 0.5*c3*y2 + 0.5*b3*x3 + 0.5*c3*y3)

    #     data = [N1*N1, 0, N1*N2, 0, N1*N3, 0,
    #             0, N1*N1, 0, N1*N2, 0, N1*N3,
    #             N2*N1, 0, N2*N2, 0, N2*N3, 0,
    #             0, N2*N1, 0, N2*N2, 0, N2*N3,
    #             N3*N1, 0, N3*N2, 0, N3*N3, 0,
    #             0, N3*N1, 0, N3*N2, 0, N3*N3]

    #     mass_matrix = self.mass_density * self.thickness * npy.array(data).reshape(6, 6)

    #     return mass_matrix.flatten()

    def strain(self):
        b_matrix = self.b_matrix()
        q = self.displacements

        strain = npy.matmul(b_matrix, q)

        return strain

    def stress(self):
        b_matrix = self.b_matrix()
        d_matrix = self.d_matrix()
        q = self.displacements

        stress = (npy.matmul(npy.matmul(d_matrix, b_matrix), q))

        return stress


class ElasticityTetrahedralElement3D(ElasticityElement, vmmesh.TetrahedralElement):
    # _standalone_in_db = False
    # _non_serializable_attributes = []
    # _non_eq_attributes = ['name']
    # _non_hash_attributes = ['name']
    # _generic_eq = True
    # def __init__(self, mesh_element: vmmesh.TetrahedralElement,
    #              elasticity_modulus, poisson_ratio,
    #              displacements = None,
    #              stress = None,
    #              strain = None,
    #              name : str = ''):
    #     self.mesh_element = mesh_element
    #     self.elasticity_modulus = elasticity_modulus
    #     self.poisson_ratio = poisson_ratio
    #     self.points = self.mesh_element.points
    #     self.b_matrix = self._b_matrix()
    #     self.d_matrix = self._d_matrix()
    #     self.displacements = displacements
    #     self.stress = stress
    #     self.strain = strain

    def __init__(self, mesh_element: vmmesh.TetrahedralElement,
                 elasticity_modulus, poisson_ratio, mass_density,
                 displacements = None,
                 stress = None,
                 strain = None,
                 name : str = ''):

        ElasticityElement.__init__(self, mesh_element,
                                       elasticity_modulus, poisson_ratio, mass_density)
        vmmesh.TetrahedralElement.__init__(self, points=mesh_element.points)

        # DessiaObject.__init__(self, name=name)

    def _b_matrix(self):

        N1, N2, N3, N4 = self.mesh_element.form_functions

        a1 = N1[1]
        b1 = N1[2]
        c1 = N1[3]
        a2 = N2[1]
        b2 = N2[2]
        c2 = N2[3]
        a3 = N3[1]
        b3 = N3[2]
        c3 = N3[3]
        a4 = N4[1]
        b4 = N4[2]
        c4 = N4[3]

        data = [a1, 0, 0, a2, 0, 0, a3, 0, 0, a4, 0, 0,
                0, b1, 0, 0, b2, 0, 0, b3, 0, 0, b4, 0,
                0, 0, c1, 0, 0, c2, 0, 0, c3, 0, 0, c4,
                b1, a1, 0, b2, a2, 0, b3, a3, 0, b4, a4, 0,
                0, c1, b1, 0, c2, b2, 0, c3, b3, 0, c4, b4,
                c1, 0, a1, c2, 0, a2, c3, 0, a3, c4, 0, a4]

        b_matrix = (1/(6*self.mesh_element.volume)) * npy.array(data).reshape(6, 12)

        return b_matrix

    def _d_matrix_plane_strain(self):
        return self._d_matrix_()
    def _d_matrix_plane_stress(self):
        return self._d_matrix_()

    def _d_matrix_(self):
        elasticity_modulus = self.elasticity_modulus
        poisson_ratio = self.poisson_ratio
        coeff_a = 1 - poisson_ratio
        coeff_b = (1 - 2 * poisson_ratio)/2

        data = [coeff_a, poisson_ratio, poisson_ratio, 0, 0, 0,
                poisson_ratio, coeff_a, poisson_ratio, 0, 0, 0,
                poisson_ratio, poisson_ratio, coeff_a, 0, 0, 0,
                0, 0, 0, coeff_b, 0, 0,
                0, 0, 0, 0, coeff_b, 0,
                0, 0, 0, 0, 0, coeff_b]

        coeff = (elasticity_modulus/((1 + poisson_ratio) * (1 - 2*poisson_ratio)))
        d_matrix = coeff * npy.array(data).reshape(6,6)

        return d_matrix

    @property
    def dimension(self):
        return 3

    def elementary_matrix(self, plane_strain: bool, plane_stress:bool):

        b_matrix = self.b_matrix
        d_matrix = self.d_matrix(plane_strain, plane_stress)

        stiffness_matrix = self.volume * (
            npy.matmul(npy.matmul(b_matrix.transpose(), d_matrix), b_matrix))

        return stiffness_matrix.flatten()

    def elementary_mass_matrix(self):
        data = [2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
                0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
                0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1,
                1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0,
                0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0,
                0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1,
                1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0,
                0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0,
                0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1,
                1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0,
                0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0,
                0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2]

        # jacobian = []
        # for p in self.points:
        #     jacobian.extend([1, *p])
        # det_jacobian = npy.linalg.det(npy.array(jacobian).reshape(4, 4))

        det_jacobian = self.volume * 6

        mass_matrix = det_jacobian * ((self.mass_density * self.volume) / 20) \
            * npy.array(data).reshape(12, 12)

        # mass_matrix = ((self.mass_density * self.volume) / 20) \
        #     * npy.array(data).reshape(12, 12)

        return mass_matrix.flatten()
