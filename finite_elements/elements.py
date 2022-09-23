#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing objects related to different finite elements types
"""

import volmdlr.mesh as vmmesh


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

    # def elementary_source_matrix(self, indexes):
    #     """
    #     Create the elementary source matrix of the MagneticElement2D

    #     :return: (double_integral_N1_dS, double_integral_N2_dS, double_integral_N3_dS)
    #     """

    #     x1 = self.triangular_element.points[0][0]
    #     y1 = self.triangular_element.points[0][1]
    #     x2 = self.triangular_element.points[1][0]
    #     y2 = self.triangular_element.points[1][1]
    #     x3 = self.triangular_element.points[2][0]
    #     y3 = self.triangular_element.points[2][1]

    #     det_jacobien = abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))

    #     element_form_functions = self.triangular_element.form_functions
    #     a1 = element_form_functions[0][0]
    #     b1 = element_form_functions[0][1]
    #     c1 = element_form_functions[0][2]
    #     a2 = element_form_functions[1][0]
    #     b2 = element_form_functions[1][1]
    #     c2 = element_form_functions[1][2]
    #     a3 = element_form_functions[2][0]
    #     b3 = element_form_functions[2][1]
    #     c3 = element_form_functions[2][2]

    #     double_integral_N1_dS = det_jacobien*(a1 + 0.5*b1*x2 + 0.5*c1*y2 + 0.5*b1*x3 + 0.5*c1*y3)
    #     double_integral_N2_dS = det_jacobien*(a2 + 0.5*b2*x2 + 0.5*c2*y2 + 0.5*b2*x3 + 0.5*c2*y3)
    #     double_integral_N3_dS = det_jacobien*(a3 + 0.5*b3*x2 + 0.5*c3*y2 + 0.5*b3*x3 + 0.5*c3*y3)

    #     return (double_integral_N1_dS, double_integral_N2_dS, double_integral_N3_dS)

    def element_to_node_factors(self):
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

