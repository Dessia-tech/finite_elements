#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing objects related to different finite elements types
"""

import volmdlr.mesh as vmmesh
from dessia_common import DessiaObject
import numpy as npy
import finite_elements.core


class Element2D(vmmesh.TriangularElement2D):
    """
    This class
    """

    def element_to_node_factors(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        x_1 = self.mesh_element.points[0][0]
        y_1 = self.mesh_element.points[0][1]
        x_2 = self.mesh_element.points[1][0]
        y_2 = self.mesh_element.points[1][1]
        x_3 = self.mesh_element.points[2][0]
        y_3 = self.mesh_element.points[2][1]

        det_jacobien = abs((x_2 - x_1) * (y_3 - y_1) - (x_3 - x_1) * (y_2 - y_1))

        element_form_functions = self.mesh_element.form_functions
        a_1 = element_form_functions[0][0]
        b_1 = element_form_functions[0][1]
        c_1 = element_form_functions[0][2]
        a_2 = element_form_functions[1][0]
        b_2 = element_form_functions[1][1]
        c_2 = element_form_functions[1][2]
        a_3 = element_form_functions[2][0]
        b_3 = element_form_functions[2][1]
        c_3 = element_form_functions[2][2]

        double_integral_n1_ds = det_jacobien * (a_1 + 0.5 * b_1 * x_2 + 0.5 * c_1 * y_2 +
                                                0.5 * b_1 * x_3 + 0.5 * c_1 * y_3)
        double_integral_n2_ds = det_jacobien * (a_2 + 0.5 * b_2 * x_2 + 0.5 * c_2 * y_2 +
                                                0.5 * b_2 * x_3 + 0.5 * c_2 * y_3)
        double_integral_n3_ds = det_jacobien * (a_3 + 0.5 * b_3 * x_2 + 0.5 * c_3 * y_2 +
                                                0.5 * b_3 * x_3 + 0.5 * c_3 * y_3)

        return (double_integral_n1_ds, double_integral_n2_ds, double_integral_n3_ds)


class MagneticElement2D(Element2D):
    """
    This class

    :param triangular_element: DESCRIPTION
    :type triangular_element: vmmesh.TriangularElement2D
    :param mu_total: DESCRIPTION
    :type mu_total: float
    :param name: DESCRIPTION, defaults to ''
    :type name: str, optional
    """

    # _standalone_in_db = False
    # _non_serializable_attributes = []
    # _non_eq_attributes = ['name']
    # _non_hash_attributes = ['name']
    # _generic_eq = True

    def __init__(self, triangular_element: vmmesh.TriangularElement2D,
                 mu_total: float, name: str = ''):
        self.triangular_element = triangular_element
        vmmesh.TriangularElement2D.__init__(self, points=triangular_element.points)
        self.mu_total = mu_total

        # DessiaObject.__init__(self, name=name)

    @property
    def dimension(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return 1

    def elementary_matrix(self):
        """
        Create the elementary matrix of the MagneticElement2D

        :return: data
        """

        element_form_functions = self.triangular_element.form_functions
        b_1 = element_form_functions[0][1]
        c_1 = element_form_functions[0][2]
        b_2 = element_form_functions[1][1]
        c_2 = element_form_functions[1][2]
        b_3 = element_form_functions[2][1]
        c_3 = element_form_functions[2][2]

        data = (1 / self.mu_total * (b_1**2 + c_1**2) * self.triangular_element.area,
                1 / self.mu_total * (b_1 * b_2 + c_1 * c_2) * self.triangular_element.area,
                1 / self.mu_total * (b_1 * b_3 + c_1 * c_3) * self.triangular_element.area,
                1 / self.mu_total * (b_1 * b_2 + c_1 * c_2) * self.triangular_element.area,
                1 / self.mu_total * (b_2**2 + c_2**2) * self.triangular_element.area,
                1 / self.mu_total * (b_2 * b_3 + c_2 * c_3) * self.triangular_element.area,
                1 / self.mu_total * (b_1 * b_3 + c_1 * c_3) * self.triangular_element.area,
                1 / self.mu_total * (b_2 * b_3 + c_2 * c_3) * self.triangular_element.area,
                1 / self.mu_total * (b_3**2 + c_3**2) * self.triangular_element.area)

        return data

    # def elementary_source_matrix(self, indexes):
    #     """
    #     Create the elementary source matrix of the MagneticElement2D

    #     :return: (double_integral_n1_ds, double_integral_n2_ds, double_integral_n3_ds)
    #     """

    #     x_1 = self.triangular_element.points[0][0]
    #     y_1 = self.triangular_element.points[0][1]
    #     x_2 = self.triangular_element.points[1][0]
    #     y_2 = self.triangular_element.points[1][1]
    #     x_3 = self.triangular_element.points[2][0]
    #     y_3 = self.triangular_element.points[2][1]

    #     det_jacobien = abs((x_2-x_1)*(y_3-y_1) - (x_3-x_1)*(y_2-y_1))

    #     element_form_functions = self.triangular_element.form_functions
    #     a_1 = element_form_functions[0][0]
    #     b_1 = element_form_functions[0][1]
    #     c_1 = element_form_functions[0][2]
    #     a_2 = element_form_functions[1][0]
    #     b_2 = element_form_functions[1][1]
    #     c_2 = element_form_functions[1][2]
    #     a_3 = element_form_functions[2][0]
    #     b_3 = element_form_functions[2][1]
    #     c_3 = element_form_functions[2][2]

    #     double_integral_n1_ds = det_jacobien*(a_1 + 0.5*b_1*x_2 + 0.5*c_1*y_2 + \
    #                                           0.5*b_1*x_3 + 0.5*c_1*y_3)
    #     double_integral_n2_ds = det_jacobien*(a_2 + 0.5*b_2*x_2 + 0.5*c_2*y_2 + \
    #                                           0.5*b_2*x_3 + 0.5*c_2*y_3)
    #     double_integral_n3_ds = det_jacobien*(a_3 + 0.5*b_3*x_2 + 0.5*c_3*y_2 + \
    #                                           0.5*b_3*x_3 + 0.5*c_3*y_3)

    #     return (double_integral_n1_ds, double_integral_n2_ds, double_integral_n3_ds)

    def element_to_node_factors(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        x_1 = self.triangular_element.points[0][0]
        y_1 = self.triangular_element.points[0][1]
        x_2 = self.triangular_element.points[1][0]
        y_2 = self.triangular_element.points[1][1]
        x_3 = self.triangular_element.points[2][0]
        y_3 = self.triangular_element.points[2][1]

        det_jacobien = abs((x_2 - x_1) * (y_3 - y_1) - (x_3 - x_1) * (y_2 - y_1))

        element_form_functions = self.triangular_element.form_functions
        a_1 = element_form_functions[0][0]
        b_1 = element_form_functions[0][1]
        c_1 = element_form_functions[0][2]
        a_2 = element_form_functions[1][0]
        b_2 = element_form_functions[1][1]
        c_2 = element_form_functions[1][2]
        a_3 = element_form_functions[2][0]
        b_3 = element_form_functions[2][1]
        c_3 = element_form_functions[2][2]

        double_integral_n1_ds = det_jacobien * (a_1 + 0.5 * b_1 * x_2 + 0.5 * c_1 * y_2 +
                                                0.5 * b_1 * x_3 + 0.5 * c_1 * y_3)
        double_integral_n2_ds = det_jacobien * (a_2 + 0.5 * b_2 * x_2 + 0.5 * c_2 * y_2 +
                                                0.5 * b_2 * x_3 + 0.5 * c_2 * y_3)
        double_integral_n3_ds = det_jacobien * (a_3 + 0.5 * b_3 * x_2 + 0.5 * c_3 * y_2 +
                                                0.5 * b_3 * x_3 + 0.5 * c_3 * y_3)

        return (double_integral_n1_ds, double_integral_n2_ds, double_integral_n3_ds)


class ElasticityElement(DessiaObject):
    """
    This class

    :param mesh_element: DESCRIPTION
    :type mesh_element: TYPE
    :param elasticity_modulus: DESCRIPTION
    :type elasticity_modulus: TYPE
    :param poisson_ratio: DESCRIPTION
    :type poisson_ratio: TYPE
    :param mass_density: DESCRIPTION
    :type mass_density: TYPE
    :param displacements: DESCRIPTION, defaults to None
    :type displacements: TYPE, optional
    :param stress: DESCRIPTION, defaults to None
    :type stress: TYPE, optional
    :param strain: DESCRIPTION, defaults to None
    :type strain: TYPE, optional
    :param name: DESCRIPTION, defaults to ''
    :type name: str, optional
    """

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
                 displacements=None,
                 stress=None,
                 strain=None,
                 name: str = ''):
        self.mesh_element = mesh_element
        self.elasticity_modulus = elasticity_modulus
        self.poisson_ratio = poisson_ratio
        self.mass_density = mass_density
        self.points = self.mesh_element.points
        self.b_matrix = self._b_matrix()
        self.d_matrix_plane_strain = self._d_matrix_plane_strain()
        self.d_matrix_plane_stress = self._d_matrix_plane_stress()
        self.displacements = displacements
        self.stress = stress
        self.strain = strain
        self.name = name

        DessiaObject.__init__(self, name=name)

    def d_matrix(self, plane_strain: bool, plane_stress: bool):
        """
        Defines

        :param plane_strain: DESCRIPTION
        :type plane_strain: bool
        :param plane_stress: DESCRIPTION
        :type plane_stress: bool
        :raises ValueError: DESCRIPTION

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if (plane_strain and plane_stress):
            raise ValueError('just one of plane_strain or plane_stress can be True')
        elif (not plane_strain and not plane_stress):
            raise ValueError('one of plane_strain or plane_stress must be True')
        elif plane_strain:
            return self.d_matrix_plane_strain
        elif plane_stress:
            return self.d_matrix_plane_stress

    def energy(self, plane_strain: bool, plane_stress: bool):
        """
        Defines

        :param plane_strain: DESCRIPTION
        :type plane_strain: bool
        :param plane_stress: DESCRIPTION
        :type plane_stress: bool

        :return: DESCRIPTION
        :rtype: TYPE
        """

        shape = self.dimension * len(self.mesh_element.points)
        return 0.5 * (npy.matmul(
            npy.matmul(npy.transpose(npy.array(self.displacements)),
                       self.elementary_matrix(plane_strain, plane_stress).reshape(shape, shape)),
            npy.array(self.displacements)))

    @classmethod
    def with_material_object(cls, mesh_element,
                             material: finite_elements.core.Material,
                             displacements=None,
                             stress=None,
                             strain=None,
                             name: str = ''):
        """
        Defines

        :param cls: DESCRIPTION
        :type cls: TYPE
        :param mesh_element: DESCRIPTION
        :type mesh_element: TYPE
        :param material: DESCRIPTION
        :type material: finite_elements.core.Material
        :param displacements: DESCRIPTION, defaults to None
        :type displacements: TYPE, optional
        :param stress: DESCRIPTION, defaults to None
        :type stress: TYPE, optional
        :param strain: DESCRIPTION, defaults to None
        :type strain: TYPE, optional
        :param name: DESCRIPTION, defaults to ''
        :type name: str, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return cls(mesh_element=mesh_element,
                   elasticity_modulus=material.elasticity_modulus,
                   poisson_ratio=material.poisson_ratio,
                   mass_density=material.mass_density,
                   displacements=displacements,
                   stress=stress,
                   strain=strain,
                   name=name)


class ElasticityTriangularElement2D(ElasticityElement, Element2D):
    """
    This class

    :param mesh_element: DESCRIPTION
    :type mesh_element: vmmesh.TriangularElement2D
    :param elasticity_modulus: DESCRIPTION
    :type elasticity_modulus: TYPE
    :param poisson_ratio: DESCRIPTION
    :type poisson_ratio: TYPE
    :param mass_density: DESCRIPTION
    :type mass_density: TYPE
    :param thickness: DESCRIPTION, defaults to 1.0
    :type thickness: float, optional
    :param displacements: DESCRIPTION, defaults to None
    :type displacements: TYPE, optional
    :param stress: DESCRIPTION, defaults to None
    :type stress: TYPE, optional
    :param strain: DESCRIPTION, defaults to None
    :type strain: TYPE, optional
    :param name: DESCRIPTION, defaults to ''
    :type name: str, optional
    """

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

    def __init__(self, mesh_element: vmmesh.TriangularElement2D,
                 elasticity_modulus, poisson_ratio, mass_density,
                 thickness: float = 1.0,
                 displacements=None,
                 stress=None,
                 strain=None,
                 name: str = ''):
        self.thickness = thickness

        ElasticityElement.__init__(self, mesh_element,
                                   elasticity_modulus, poisson_ratio, mass_density, name=name)
        vmmesh.TriangularElement2D.__init__(self, points=mesh_element.points, name=name)

        # DessiaObject.__init__(self, name=name)

    def _b_matrix(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        y_coor = [(self.points[i].y - self.points[j].y) for (i, j) in [(1, 2), (2, 0), (0, 1)]]
        x_coor = [(self.points[i].x - self.points[j].x) for (i, j) in [(2, 1), (0, 2), (1, 0)]]

        det_jacobian = \
            (self.points[0].x - self.points[2].x) * (self.points[1].y - self.points[2].y) \
            - (self.points[0].y - self.points[2].y) * (self.points[1].x - self.points[2].x)

        data = [y_coor[0], 0, y_coor[1], 0, y_coor[2], 0,
                0, x_coor[0], 0, x_coor[1], 0, x_coor[2],
                x_coor[0], y_coor[0], x_coor[1], y_coor[1], x_coor[2], y_coor[2]]

        b_matrix = (1 / det_jacobian) * npy.array(data).reshape(3, 6)

        return b_matrix

    def _d_matrix_plane_strain(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        a_param = (self.elasticity_modulus * self.poisson_ratio) \
            / ((1 + self.poisson_ratio) * (1 - 2 * self.poisson_ratio))
        b_param = self.elasticity_modulus / (2 * (1 + self.poisson_ratio))
        c_param = a_param + 2 * b_param
        data = [c_param, a_param, 0,
                a_param, c_param, 0,
                0, 0, b_param]

        return npy.array(data).reshape(3, 3)

    def _d_matrix_plane_stress(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        a_param = (self.elasticity_modulus * self.poisson_ratio) \
            / (1 - self.poisson_ratio**2)
        b_param = self.elasticity_modulus / (2 * (1 + self.poisson_ratio))
        c_param = a_param + 2 * b_param

        data = [c_param, a_param, 0,
                a_param, c_param, 0,
                0, 0, b_param]

        return npy.array(data).reshape(3, 3)

    @property
    def dimension(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return 2

    def elementary_matrix(self, plane_strain: bool, plane_stress: bool):
        """
        Defines

        :param plane_strain: DESCRIPTION
        :type plane_strain: bool
        :param plane_stress: DESCRIPTION
        :type plane_stress: bool

        :return: DESCRIPTION
        :rtype: TYPE
        """

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

    @classmethod
    def from_element(cls, mesh_element, elasticity_element):
        """
        Defines

        :param cls: DESCRIPTION
        :type cls: TYPE
        :param mesh_element: DESCRIPTION
        :type mesh_element: TYPE
        :param elasticity_element: DESCRIPTION
        :type elasticity_element: TYPE

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return cls(mesh_element, elasticity_element.elasticity_modulus,
                   elasticity_element.poisson_ratio, elasticity_element.mass_density,
                   elasticity_element.thickness)

    def strain(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        b_matrix = self.b_matrix()
        q_vector = self.displacements

        strain = npy.matmul(b_matrix, q_vector)

        return strain

    def stress(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        b_matrix = self.b_matrix()
        d_matrix = self.d_matrix()
        q_vector = self.displacements

        stress = (npy.matmul(npy.matmul(d_matrix, b_matrix), q_vector))

        return stress


class ElasticityTetrahedralElement3D(ElasticityElement, vmmesh.TetrahedralElement):
    """
    This class

    :param mesh_element: DESCRIPTION
    :type mesh_element: vmmesh.TetrahedralElement
    :param elasticity_modulus: DESCRIPTION
    :type elasticity_modulus: TYPE
    :param poisson_ratio: DESCRIPTION
    :type poisson_ratio: TYPE
    :param mass_density: DESCRIPTION
    :type mass_density: TYPE
    :param displacements: DESCRIPTION, defaults to None
    :type displacements: TYPE, optional
    :param stress: DESCRIPTION, defaults to None
    :type stress: TYPE, optional
    :param strain: DESCRIPTION, defaults to None
    :type strain: TYPE, optional
    :param name: DESCRIPTION, defaults to ''
    :type name: str, optional
    """

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
                 displacements=None,
                 stress=None,
                 strain=None,
                 name: str = ''):

        ElasticityElement.__init__(self, mesh_element,
                                   elasticity_modulus, poisson_ratio, mass_density, name=name)
        vmmesh.TetrahedralElement.__init__(self, points=mesh_element.points, name=name)

        # DessiaObject.__init__(self, name=name)

    def _b_matrix(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        form_fct_1, form_fct_2, form_fct_3, form_fct_4 = self.mesh_element.form_functions

        a_1 = form_fct_1[1]
        b_1 = form_fct_1[2]
        c_1 = form_fct_1[3]
        a_2 = form_fct_2[1]
        b_2 = form_fct_2[2]
        c_2 = form_fct_2[3]
        a_3 = form_fct_3[1]
        b_3 = form_fct_3[2]
        c_3 = form_fct_3[3]
        a_4 = form_fct_4[1]
        b_4 = form_fct_4[2]
        c_4 = form_fct_4[3]

        data = [a_1, 0, 0, a_2, 0, 0, a_3, 0, 0, a_4, 0, 0,
                0, b_1, 0, 0, b_2, 0, 0, b_3, 0, 0, b_4, 0,
                0, 0, c_1, 0, 0, c_2, 0, 0, c_3, 0, 0, c_4,
                b_1, a_1, 0, b_2, a_2, 0, b_3, a_3, 0, b_4, a_4, 0,
                0, c_1, b_1, 0, c_2, b_2, 0, c_3, b_3, 0, c_4, b_4,
                c_1, 0, a_1, c_2, 0, a_2, c_3, 0, a_3, c_4, 0, a_4]

        b_matrix = (1 / (6 * self.mesh_element.volume)) * npy.array(data).reshape(6, 12)

        return b_matrix

    def _d_matrix_plane_strain(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return self._d_matrix_()

    def _d_matrix_plane_stress(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return self._d_matrix_()

    def _d_matrix_(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        elasticity_modulus = self.elasticity_modulus
        poisson_ratio = self.poisson_ratio
        coeff_a = 1 - poisson_ratio
        coeff_b = (1 - 2 * poisson_ratio) / 2

        data = [coeff_a, poisson_ratio, poisson_ratio, 0, 0, 0,
                poisson_ratio, coeff_a, poisson_ratio, 0, 0, 0,
                poisson_ratio, poisson_ratio, coeff_a, 0, 0, 0,
                0, 0, 0, coeff_b, 0, 0,
                0, 0, 0, 0, coeff_b, 0,
                0, 0, 0, 0, 0, coeff_b]

        coeff = (elasticity_modulus / ((1 + poisson_ratio) * (1 - 2 * poisson_ratio)))
        d_matrix = coeff * npy.array(data).reshape(6, 6)

        return d_matrix

    @property
    def dimension(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return 3

    def elementary_matrix(self, plane_strain: bool, plane_stress: bool):
        """
        Defines

        :param plane_strain: DESCRIPTION
        :type plane_strain: bool
        :param plane_stress: DESCRIPTION
        :type plane_stress: bool

        :return: DESCRIPTION
        :rtype: TYPE
        """

        b_matrix = self.b_matrix
        d_matrix = self.d_matrix(plane_strain, plane_stress)

        stiffness_matrix = self.volume * (
            npy.matmul(npy.matmul(b_matrix.transpose(), d_matrix), b_matrix))

        return stiffness_matrix.flatten()

    @classmethod
    def from_element(cls, mesh_element, elasticity_element):
        """
        Defines

        :param cls: DESCRIPTION
        :type cls: TYPE
        :param mesh_element: DESCRIPTION
        :type mesh_element: TYPE
        :param elasticity_element: DESCRIPTION
        :type elasticity_element: TYPE

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return cls(mesh_element, elasticity_element.elasticity_modulus,
                   elasticity_element.poisson_ratio, elasticity_element.mass_density)
