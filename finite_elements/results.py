#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing objects related to finite elements analysis results
"""

from typing import List  # Tuple, TypeVar
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner
# from matplotlib import cm
import numpy as npy
import volmdlr as vm
import volmdlr.mesh as vmmesh
from dessia_common import DessiaObject
from finite_elements.core import MU, blue_red
import finite_elements.core


class Result(DessiaObject):
    """
    This class

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

    @property
    def dimension(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return self.mesh.elements_groups[0].elements[0].dimension


class MagneticResults(Result):
    """
    This class

    :param mesh: DESCRIPTION
    :type mesh: vmmesh.Mesh
    :param result_vector: DESCRIPTION
    :type result_vector: List[float]
    """

    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, mesh: vmmesh.Mesh, result_vector: List[float]):
        self.mesh = mesh
        self.result_vector = result_vector

        self.magnetic_field_per_element = self._magnetic_field_per_element()
        self.magnetic_field_norm = self._magnetic_field_norm()

        Result.__init__(self, mesh, result_vector)

    def brbtetha(self, gap_elements_group):
        """
        Defines

        :param gap_elements_group: DESCRIPTION
        :type gap_elements_group: TYPE

        :return: DESCRIPTION
        :rtype: TYPE
        """

        element_to_magnetic_field = self.magnetic_field_per_element
        all_br_btetha = []
        for element in gap_elements_group.elements:
            vector_b = element_to_magnetic_field[element]

            element_center = element.center
            e_r = vm.Vector2D(element_center.vector)
            e_r.Normalize()
            e_teta = vm.Vector2D((-e_r[1], e_r[0]))
            b_r = vector_b.Dot(e_r)
            b_teta = vector_b.Dot(e_teta)

            all_br_btetha.append(b_r * b_teta)
        return all_br_btetha

    def _magnetic_field_norm(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        element_to_magnetic_field = self.magnetic_field_per_element
        bs_param = [b_param.norm() for b_param in list(element_to_magnetic_field.values())]

        return bs_param

    def _magnetic_field_per_element(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        element_to_magnetic_field = {}
        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:
                element_form_functions = element.form_functions
                indexes = [self.mesh.node_to_index[element.points[0]],
                           self.mesh.node_to_index[element.points[1]],
                           self.mesh.node_to_index[element.points[2]]]
                b1_param = element_form_functions[0][1]
                c1_param = element_form_functions[0][2]
                b2_param = element_form_functions[1][1]
                c2_param = element_form_functions[1][2]
                b3_param = element_form_functions[2][1]
                c3_param = element_form_functions[2][2]
                b_x = float(c1_param *
                            self.result_vector[indexes[0]] +
                            c2_param *
                            self.result_vector[indexes[1]] +
                            c3_param *
                            self.result_vector[indexes[2]])
                b_y = float(-b1_param * self.result_vector[indexes[0]] - b2_param *
                            self.result_vector[indexes[1]] - b3_param
                            * self.result_vector[indexes[2]])
                element_to_magnetic_field[element] = vm.Vector2D(b_x, b_y)
        return element_to_magnetic_field

    def maxwell_stress_tensor(self, element):
        """
        Computes the Maxwell stress tensor for one element.
            Returns a list made of sigma_r_r, sigma_r_theta and sigma_theta_theta,
            since sigma_r_theta = sigma_theta_r

        :param element: The element on which the tensor is computed.
        :type element: a TriangularElement object

        :return: DESCRIPTION
        :rtype: TYPE
        """

        element_to_magnetic_field = self.magnetic_field_per_element
        vector_b = element_to_magnetic_field[element]
        element_center = element.center
        e_r = vm.Vector2D(element_center.vector)
        e_r.Normalize()
        e_teta = vm.Vector2D((-e_r[1], e_r[0]))
        b_r = vector_b.Dot(e_r)
        b_teta = vector_b.Dot(e_teta)

        sigma_rr = 1 / MU * b_r**2 - 1 / (2 * MU) * vector_b.Norm()**2
        sigma_rteta = 1 / MU * b_r * b_teta
        sigma_tetateta = 1 / MU * b_teta**2 - 1 / (2 * MU) * vector_b.Norm()**2
        sigma_rr_rteta_tetateta = [sigma_rr, sigma_rteta, sigma_tetateta]
        return sigma_rr_rteta_tetateta

    def torque(self, air_gap_elements_group_name, length_motor, radius_stator, radius_rotor):
        # nb_notches):
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

        :return: DESCRIPTION
        :rtype: TYPE
        """

        element_to_magnetic_field = self.magnetic_field_per_element

        for elements_group in self.mesh.elements_groups:
            if elements_group.name == air_gap_elements_group_name:
                gap_elements_group = elements_group
                break

        somme = 0
        i = 0

        # r = (radius_stator - radius_rotor)/2 + radius_rotor

        # fig, ax = plt.subplots()

        for element in gap_elements_group.elements:
            vector_b = element_to_magnetic_field[element]

            element_center = element.center
            e_r = vm.Vector2D(element_center.vector)
            e_r.Normalize()
            e_teta = vm.Vector2D((-e_r[1], e_r[0]))
            b_r = vector_b.Dot(e_r)
            b_teta = vector_b.Dot(e_teta)
            r_br_bteta = element_center.Norm() * b_r * b_teta
            # r_br_bteta = r * b_r * b_teta
            ds_param = element.area

            # e_r.plot(ax=ax, origin=element_center, amplitude=0.005, color='b')
            # e_teta.plot(ax=ax, origin=element_center, amplitude=0.005, color='g')
            # vector_b.plot(ax=ax, origin=element_center, amplitude=0.005, color='r')

            somme += r_br_bteta * ds_param
            i += 1

        # print('nb elements in airgap', i)
        t_param = length_motor / (MU * (radius_stator - radius_rotor)) * somme

        return t_param

    def plot_brbtetha(self, ax=None, air_gap_elements_group_name='Gap ring'):
        """
        Plots

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional
        :param air_gap_elements_group_name: DESCRIPTION, defaults to 'Gap ring'
        :type air_gap_elements_group_name: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

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

        all_br_btetha = self.brbtetha(gap_elements_group)

        # color_map = ((0,0,1), (1,0,0))
        jet = plt.get_cmap('jet')
        # bs_param = [B.Norm() for B in list(element_to_magnetic_field.values())]

        br_btetha_max, br_btetha_min = max(all_br_btetha), min(all_br_btetha)

        b_to_color = {}
        all_colors = []
        for b_param in all_br_btetha:
            if b_param > br_btetha_max:
                x_param = 1
            else:
                x_param = (b_param - br_btetha_min) / (br_btetha_max - br_btetha_min)
            # color = (color_map[0][0]-(color_map[0][0]-color_map[1][0])*x,
            #          color_map[0][1]-(color_map[0][1]-color_map[1][1])*x,
            #          color_map[0][2]-(color_map[0][2]-color_map[1][2])*x)
            color = jet(int(x_param * 256))[:3]

            b_to_color[b_param] = color
            all_colors.append(color)
        # print(b_to_color)
        for i, element in enumerate(gap_elements_group.elements):
            # element.plot(ax=ax, color=b_to_color[element_to_magnetic_field[element].Norm()],
            #              fill=True)
            element.plot(ax=ax, color=all_colors[i], fill=True)

        norm = mpl.colors.Normalize(vmin=br_btetha_min, vmax=br_btetha_max)
        scalar_mappable = plt.cm.ScalarMappable(cmap=jet, norm=norm)

        scalar_mappable.set_array([])
        cbar = fig.colorbar(scalar_mappable, ticks=npy.linspace(br_btetha_min, br_btetha_max, 10))
        # cbar = fig.colorbar(sm, ticks=npy.linspace(-0.9, 0.8, 10))
        cbar.set_label('Br*Btetha')

        return ax

    def plot_magnetic_field(self, ax=None, bmax=None):
        """
        Plots the mesh with colored triangular elements representing the \
        intensity of the induction field inside the Machine.

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional
        :param bmax: DESCRIPTION, defaults to None
        :type bmax: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()
        element_to_magnetic_field = self.magnetic_field_per_element
        bs_param = self.magnetic_field_norm

        b_max, b_min = finite_elements.core.get_bmin_bmax(bs_param,
                                                          param_bmin=None,
                                                          param_bmax=bmax)
        b_to_color = finite_elements.core.get_colors(bs_param, param_bmax=b_max, param_bmin=b_min)

        for group in self.mesh.elements_groups:
            for element in group.elements:
                element.plot(ax=ax, color=b_to_color[element_to_magnetic_field[element].norm()],
                             fill=True)

        norm = mpl.colors.Normalize(vmin=b_min, vmax=b_max)
        scalar_mappable = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)
        scalar_mappable.set_array([])
        cbar = fig.colorbar(scalar_mappable, ticks=npy.linspace(b_min, b_max, 10))
        cbar.set_label('Magnetic Field in Tesla')

        return ax

    def plot_magnetic_field_contour(self, ax=None):
        """
        Plots

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if ax is None:
            _, ax = plt.subplots()
        else:
            _ = plt.gcf()
        ax.set_aspect('equal')

        element_to_magnetic_field = self.magnetic_field_per_element

        x_coor = []
        y_coor = []
        z_param = []
        for group in self.mesh.elements_groups:
            for element in group.elements:
                x_center, y_center = element.center
                x_coor.append(x_center)
                y_coor.append(y_center)
                z_param.append(element_to_magnetic_field[element].norm())

        tri = Triangulation(x_coor, y_coor)

        # -----------------------------------------------------------------------------
        # Improving the triangulation before high-res plots: removing flat triangles
        # -----------------------------------------------------------------------------
        # masking badly shaped triangles at the border of the triangular mesh.
        min_circle_ratio = -1
        mask = TriAnalyzer(tri).get_flat_tri_mask(min_circle_ratio)
        tri.set_mask(mask)

        # refining the data
        refiner = UniformTriRefiner(tri)
        subdiv = 3
        tri_refi, z_test_refi = refiner.refine_field(z_param, subdiv=subdiv)

        levels = npy.arange(0., 1., 0.05)
        # cmap = cm.get_cmap(name='Blues', lut=None)
        ax.tricontour(tri_refi, z_test_refi, levels=levels,  # cmap=cmap,
                      linewidths=[2.0, 0.5, 1.0, 0.5])
        # ax.triplot(tri_refi, color='0.97')
        # ax.triplot(tri, color='0.7')
        # ax.tricontour(x, y, Z)
        return ax

    def plot_magnetic_field_vectors(self, ax=None, amplitude=0.005,
                                    bmax=None, bmin=None):
        """
        Plots the mesh with a field of vectors representing the induction field \
        and its direction inside the Machine.

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional
        :param amplitude: DESCRIPTION, defaults to 0.005
        :type amplitude: TYPE, optional
        :param bmax: DESCRIPTION, defaults to None
        :type bmax: TYPE, optional
        :param bmin: DESCRIPTION, defaults to None
        :type bmin: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()

        element_to_magnetic_field = self.magnetic_field_per_element
        bs_param = self.magnetic_field_norm

        b_max, b_min = finite_elements.core.get_bmin_bmax(param_bs=bs_param, param_bmax=bmax,
                                                          param_bmin=bmin)
        b_to_color = finite_elements.core.get_colors(bs_param, param_bmax=b_max, param_bmin=b_min)

        for element, b_param in element_to_magnetic_field.items():
            b_param.plot(amplitude=amplitude, origin=element.center, ax=ax,
                         color=b_to_color[b_param.norm()], normalize=True)

        norm = mpl.colors.Normalize(vmin=b_min, vmax=b_max)
        scalar_mappable = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)
        scalar_mappable.set_array([])
        cbar = fig.colorbar(scalar_mappable, ticks=npy.linspace(b_min, b_max, 10))
        cbar.set_label('Magnetic Field in Tesla')

        return ax

    def plot_potential_vector(self, ax=None):
        """
        Plots the mesh with colored triangular elements representing the \
        intensity of the magnetic potential inside the Machine.

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        triang = finite_elements.core.get_triangulation(self.mesh)
        # z_param = npy.asarray([p for p in self.result_vector[:len(self.mesh.nodes)]])  # p[0]
        z_param = npy.asarray(list(self.result_vector[:len(self.mesh.nodes)]))  # unnecessary-comprehension
        z_min, z_max = min(z_param), max(z_param)

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()
        ax.tricontourf(triang, z_param, cmap=blue_red)
        ax.triplot(triang, 'k-')
        ax.set_title('Triangular grid')

        norm = mpl.colors.Normalize(vmin=z_min, vmax=z_max)
        scalar_mappable = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)
        scalar_mappable.set_array([])
        cbar = fig.colorbar(scalar_mappable, ticks=npy.linspace(z_min, z_max, 10))
        cbar.set_label('Potential Vector')

        return ax


class ElasticityResults(Result):
    """
    This class

    :param mesh: DESCRIPTION
    :type mesh: vmmesh.Mesh
    :param result_vector: DESCRIPTION
    :type result_vector: List[float]
    :param plane_strain: DESCRIPTION
    :type plane_strain: bool
    :param plane_stress: DESCRIPTION
    :type plane_stress: bool
    """

    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, mesh: vmmesh.Mesh, result_vector: List[float],
                 plane_strain: bool,
                 plane_stress: bool):
        self.mesh = mesh
        self.result_vector = result_vector
        self.plane_strain = plane_strain
        self.plane_stress = plane_stress

        self.displacement_vectors_per_node = self._displacement_vectors_per_node()
        self.displacements_per_element = self._displacements_per_element()
        self.energy_per_element = self._energy_per_element()
        self.energy = self._energy()
        self.strain, self.stress = self._strain_stress_per_element()
        self.deformed_nodes = self._deformed_nodes()
        self.deformed_mesh = self._deformed_mesh()

        Result.__init__(self, mesh, result_vector)

    def _displacements_per_element(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        positions = finite_elements.core.global_matrix_positions(dimension=self.dimension,
                                                                 nodes_number=len(self.mesh.nodes))
        q_vector = self.result_vector

        displacements_per_element = {}
        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:
                displacements = []
                indexes = [self.mesh.node_to_index[point] for point in element.points]
                for index in indexes:
                    for i in range(self.dimension):
                        displacements.append(q_vector[positions[(index, i + 1)]])

                displacements_per_element[element] = displacements
                element.displacements = displacements

        return displacements_per_element

    def _displacement_vectors_per_node(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        nodes_number = len(self.mesh.nodes)
        positions = finite_elements.core.global_matrix_positions(dimension=self.dimension,
                                                                 nodes_number=nodes_number)
        displacement_field_vectors = []
        q_vector = self.result_vector

        for node in range(0, nodes_number):
            displacement = []
            for i in range(self.dimension):
                displacement.append(q_vector[positions[(node, i + 1)]])

            displacement_field_vectors.append(
                getattr(vm, f'Vector{self.__class__.__name__[-2::]}')(*displacement))
            # displacement_field_vectors.append(vm.Vector2D(*displacement))

        return displacement_field_vectors

    def _strain_stress_per_element(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        element_to_strain, element_to_stress = {}, {}
        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:
                element_to_strain[element] = (npy.matmul(element.b_matrix, element.displacements))
                element.strain = element_to_strain[element]
                element_to_stress[element] = (npy.matmul(
                    npy.matmul(element.d_matrix(plane_strain=self.plane_strain,
                                                plane_stress=self.plane_stress), element.b_matrix),
                    element.displacements))
                element.stress = element_to_stress[element]

        return element_to_strain, element_to_stress

    def _deformed_mesh(self, amplitude=1):
        """
        Defines

        :param amplitude: DESCRIPTION, defaults to 1
        :type amplitude: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if amplitude == 1:
            deformed_nodes = self.deformed_nodes
        else:
            deformed_nodes = self._deformed_nodes(amplitude=amplitude)

        group_elasticity_elments = []
        for elements_group in self.mesh.elements_groups:
            elasticity_elments = []
            for element in elements_group.elements:

                indexes = [self.mesh.node_to_index[point] for point in element.points]

                points = [deformed_nodes[index] for index in indexes]

                mesh_element = element.mesh_element.__class__(points)

                elasticity_elments.append(element.__class__.from_element(mesh_element,
                                                                         element))

            group_elasticity_elments.append(vmmesh.ElementsGroup(elasticity_elments, ''))

        mesh = vmmesh.Mesh(group_elasticity_elments)
        mesh.nodes = deformed_nodes  # Keep self.mesh order
        mesh.node_to_index = {mesh.nodes[i]: i for i in range(len(mesh.nodes))}

        return mesh

    def _deformed_nodes(self, amplitude=1):
        """
        Defines

        :param amplitude: DESCRIPTION, defaults to 1
        :type amplitude: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        displacement_field_vectors = self.displacement_vectors_per_node
        deformed_nodes = []
        for i, node in enumerate(self.mesh.nodes):
            obj = getattr(vmmesh, f'Node{self.__class__.__name__[-2::]}')
            deformed_nodes.append(
                getattr(obj, 'from_point')(node + displacement_field_vectors[i] * amplitude))

        return deformed_nodes

    def _energy(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        # # shape = len(self.mesh.nodes) * self.dimension
        # # K = self.create_matrix() (!)
        # # displacements = self.result_vector[0:shape]

        # # return 0.5 * (npy.matmul(npy.matmul(npy.transpose(npy.array(displacements)),
        # #                          K),
        # #               npy.array(displacements)))

        # energy = 0
        # for group in self.mesh.elements_groups:
        #     for element in group.elements:
        #         energy += element.energy(self.plane_strain, self.plane_stress)

        # return energy

        # return sum([value for value in self.energy_per_element.values()]) #unnecessary-comprehension
        return sum(list(self.energy_per_element.values()))

    def _energy_per_element(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        energy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                energy[element] = (element.energy(self.plane_strain, self.plane_stress))
        return energy

    def displacement_per_node_x(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return [displacement[0] for displacement in self.displacement_vectors_per_node]

    def displacement_per_node_y(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return [displacement[1] for displacement in self.displacement_vectors_per_node]

    def update_vtk_with_results(self, input_file_name, output_file_name):
        """
        Defines

        :param input_file_name: DESCRIPTION
        :type input_file_name: TYPE
        :param output_file_name: DESCRIPTION
        :type output_file_name: TYPE

        :return: DESCRIPTION
        :rtype: TYPE
        """

        with open(input_file_name) as f_in:
            with open(output_file_name, "w") as f_out:
                for line in f_in:
                    f_out.write(line)
        f_out.close()
        f_in.close()

        lines = ['POINT_DATA ' + str(len(self.mesh.nodes))]
        lines.append('SCALARS ' + 'Displacement_Magnitude float 1')
        lines.append('LOOKUP_TABLE default')
        for displacement in self.displacement_vectors_per_node:
            lines.append(str(displacement.norm()))

        lines.append('VECTORS Displacement_Vectors float')
        for displacement in self.displacement_vectors_per_node:
            line = ''
            for i in range(len([*displacement])):
                line += str(displacement[i]) + ' '
            if i == 1:
                line += '0'
            lines.append(line)
            # lines.append(str(displacement.x)+' '+str(displacement.y)+' '+str(displacement.z))

        # lines.append('CELL_DATA ' + str(104))
        # lines.append('SCALARS cell_scalars int 1')
        # lines.append('LOOKUP_TABLE default')
        # for i in range(104):
        #     lines.append(str(i))
        # lines.append('SCALARS axial_strain_x float 1')
        # lines.append('LOOKUP_TABLE default')

        # axial_strain_x = self.axial_strain_x()
        # for _, value in axial_strain_x.items():
        #     lines.append(str(value))

        with open(output_file_name, "a+") as f_out:
            # f_out.seek(0)
            # f.write('\n')
            for line in lines:
                f_out.write(line)
                f_out.write('\n')
        f_out.close()


class ElasticityResults2D(ElasticityResults):
    """
    This class
    """

    # _standalone_in_db = True
    # _non_serializable_attributes = []
    # _non_eq_attributes = ['name']
    # _non_hash_attributes = ['name']
    # _generic_eq = True
    # def __init__(self, mesh: vmmesh.Mesh, result_vector: List[float]):
    #     self.mesh = mesh
    #     self.result_vector = result_vector

    #     self.displacement_vectors_per_node = self._displacement_vectors_per_node()
    #     self.displacements_per_element = self._displacements_per_element()
    #     self.strain, self.stress = self._strain_stress_per_element()
    #     self.deformed_nodes = self._deformed_nodes()
    #     self.deformed_mesh = self._deformed_mesh()

    #     Result.__init__(self, mesh, result_vector)

    def axial_strain_x(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        strain = self.strain
        axial_strain_x = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_strain_x[element] = strain[element][0]

        return axial_strain_x

    def axial_strain_y(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        strain = self.strain
        axial_strain_y = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_strain_y[element] = strain[element][1]

        return axial_strain_y

    def axial_stress_x(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        stress = self.stress
        axial_stress_x = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_stress_x[element] = stress[element][0]

        return axial_stress_x

    def axial_stress_y(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        stress = self.stress
        axial_stress_y = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_stress_y[element] = stress[element][1]

        return axial_stress_y

    def plot_axial_strain_x(self, ax=None, fig=None):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return self.plot_constraints(constraint_name='axial_strain_x', ax=ax, fig=fig)

    def plot_axial_strain_y(self, ax=None, fig=None):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return self.plot_constraints(constraint_name='axial_strain_y', ax=ax, fig=fig)

    def plot_axial_stress_x(self, ax=None, fig=None):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return self.plot_constraints(constraint_name='axial_stress_x', ax=ax, fig=fig)

    def plot_axial_stress_y(self, ax=None, fig=None):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return self.plot_constraints(constraint_name='axial_stress_y', ax=ax, fig=fig)

    def plot_constraints(self, constraint_name: str, ax=None, fig=None):
        """
        Defines

        :param constraint_name: DESCRIPTION
        :type constraint_name: str
        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional
        :param fig: DESCRIPTION, defaults to None
        :type fig: TYPE, optional
        :raises NotImplementedError: DESCRIPTION

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if ax is None:
            fig, ax = plt.subplots()

        if hasattr(self, constraint_name):
            result = getattr(self, constraint_name)()
            # result_values = [value for value in result.values()] #unnecessary-comprehension
            result_values = list(result.values())
        else:
            raise NotImplementedError(
                f'Class {self.__class__.__name__} does not implement {constraint_name}')

        deformed_mesh = self.deformed_mesh
        b_max, b_min = finite_elements.core.get_bmin_bmax(result_values, param_bmax=None,
                                                          param_bmin=None)
        b_to_color = finite_elements.core.get_colors(result_values, param_bmax=b_max,
                                                     param_bmin=b_min)
        for g_index, group in enumerate(deformed_mesh.elements_groups):
            for e_index, element in enumerate(group.elements):
                element.plot(ax=ax,
                             color=b_to_color[result[
                                 self.mesh.elements_groups[g_index].elements[e_index]]],
                             fill=True)

        norm = mpl.colors.Normalize(vmin=b_min, vmax=b_max)
        scalar_mappable = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)
        scalar_mappable.set_array([])
        cbar = fig.colorbar(scalar_mappable, ticks=npy.linspace(b_min, b_max, 10), ax=ax)
        cbar.set_label(constraint_name)

        ax.set_xlabel('x')
        ax.set_ylabel('y')

        return ax

    def plot_deformed_mesh(self, ax=None, amplitude=1):
        """
        Defines

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional
        :param amplitude: DESCRIPTION, defaults to 1
        :type amplitude: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if ax is None:
            _, ax = plt.subplots()
            ax.set_aspect('equal')
        if amplitude != 1:
            self._deformed_mesh(amplitude=amplitude).plot(ax=ax)
        else:
            self.deformed_mesh.plot(ax=ax)
        # self.mesh.plot(ax)

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        if self.plane_stress:
            comment = '"Plane Stress"'
        else:
            comment = '"Plane Strain"'

        ax.set_title('Deformed Mesh ' + comment)

        return ax

    def plot_displacements(self, displacement_name: str, ax=None, fig=None, amplitude=1):
        """
        Defines

        :param displacement_name: DESCRIPTION
        :type displacement_name: str
        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional
        :param fig: DESCRIPTION, defaults to None
        :type fig: TYPE, optional
        :param amplitude: DESCRIPTION, defaults to 1
        :type amplitude: TYPE, optional
        :raises NotImplementedError: DESCRIPTION

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if hasattr(self, displacement_name):
            x = getattr(self, displacement_name)()
        elif displacement_name == 'displacement_per_node_xy':
            x = [displacement.norm() for displacement in self.displacement_vectors_per_node]
            displacement_name = 'displacement_magnitude_norm'
        else:
            raise NotImplementedError(
                f'Class {self.__class__.__name__} does not implement {displacement_name}')

        if amplitude != 1:
            mesh_fe = self._deformed_mesh(amplitude=amplitude)
        else:
            mesh_fe = self.deformed_mesh

        triang = finite_elements.core.get_triangulation(mesh_fe)  # self.mesh
        x_min, x_max = min(x), max(x)

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()

        ax.set_aspect('equal')
        ax.tricontourf(triang, x, cmap=blue_red)
        ax.triplot(triang, 'k-')
        # ax.set_title('Triangular grid')

        norm = mpl.colors.Normalize(vmin=x_min, vmax=x_max)
        scalar_mappable = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)
        scalar_mappable.set_array([])
        cbar = fig.colorbar(scalar_mappable, ticks=npy.linspace(x_min, x_max, 10))
        cbar.set_label(displacement_name)

        ax.set_xlabel('x')
        ax.set_ylabel('y')

        if self.plane_stress:
            comment = '"Plane Stress"'
        else:
            comment = '"Plane Strain"'
        ax.set_title(comment)

        return ax

    def plot_displacement_per_node_x(self, ax=None, amplitude=1):
        """
        Defines

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional
        :param amplitude: DESCRIPTION, defaults to 1
        :type amplitude: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return self.plot_displacements(displacement_name='displacement_per_node_x', ax=ax,
                                       amplitude=amplitude)

    def plot_displacement_per_node_xy(self, ax=None, amplitude=1):
        """
        Defines

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional
        :param amplitude: DESCRIPTION, defaults to 1
        :type amplitude: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return self.plot_displacements(displacement_name='displacement_per_node_xy', ax=ax,
                                       amplitude=amplitude)

    def plot_displacement_per_node_y(self, ax=None, amplitude=1):
        """
        Defines

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional
        :param amplitude: DESCRIPTION, defaults to 1
        :type amplitude: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return self.plot_displacements(displacement_name='displacement_per_node_y', ax=ax,
                                       amplitude=amplitude)

    def plot_displacement_vectors_per_node(self, ax=None, amplitude=0.05):
        """
        Defines

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional
        :param amplitude: DESCRIPTION, defaults to 0.05
        :type amplitude: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if ax is None:
            _, ax = plt.subplots()

        ax.set_aspect('equal')
        self.mesh.plot(ax)
        displacement_field_vectors = self.displacement_vectors_per_node
        for i, vector in enumerate(displacement_field_vectors):
            vector.plot(amplitude=amplitude, origin=self.mesh.nodes[i], ax=ax, normalize=True)

        ax.set_xlabel('x')
        ax.set_ylabel('y')

        if self.plane_stress:
            comment = 'Displacement Vectors "Plane Stress"'
        else:
            comment = 'Displacement Vectors "Plane Strain"'
        ax.set_title(comment)

    def plot_energy(self, ax=None, fig=None, amplitude=1):
        """
        Defines

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional
        :param fig: DESCRIPTION, defaults to None
        :type fig: TYPE, optional
        :param amplitude: DESCRIPTION, defaults to 1
        :type amplitude: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        result = self.energy_per_element
        # result_values = [value for value in self.energy_per_element.values()] #unnecessary-comprehension
        result_values = list(self.energy_per_element.values())

        if amplitude != 1:
            deformed_mesh = self._deformed_mesh(amplitude=amplitude)
        else:
            deformed_mesh = self.deformed_mesh

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()

        b_max, b_min = finite_elements.core.get_bmin_bmax(result_values, param_bmin=None,
                                                          param_bmax=None)
        b_to_color = finite_elements.core.get_colors(result_values, param_bmax=b_max,
                                                     param_bmin=b_min)
        for g_index, group in enumerate(deformed_mesh.elements_groups):
            for e_index, element in enumerate(group.elements):
                element.plot(ax=ax,
                             color=b_to_color[result[
                                 self.mesh.elements_groups[g_index].elements[e_index]]],
                             fill=True)

        norm = mpl.colors.Normalize(vmin=b_min, vmax=b_max)
        scalar_mappable = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)
        scalar_mappable.set_array([])
        cbar = fig.colorbar(scalar_mappable, ticks=npy.linspace(b_min, b_max, 10), ax=ax)
        cbar.set_label('Energy')
        ax.set_aspect('equal')

        ax.set_xlabel('x')
        ax.set_ylabel('y')

        if self.plane_stress:
            comment = '"Plane Stress"'
        else:
            comment = '"Plane Strain"'
        ax.set_title(comment)

        return ax

    def plot_shear_strain_xy(self, ax=None, fig=None):
        """
        Plots

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional
        :param fig: DESCRIPTION, defaults to None
        :type fig: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return self.plot_constraints(constraint_name='shear_strain_xy', ax=ax, fig=fig)

    def plot_shear_stress_xy(self, ax=None, fig=None):
        """
        Plots

        :param ax: DESCRIPTION, defaults to None
        :type ax: TYPE, optional
        :param fig: DESCRIPTION, defaults to None
        :type fig: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return self.plot_constraints(constraint_name='shear_stress_xy', ax=ax, fig=fig)

    def plot_strain(self, axs=None, fig=None, row=1):
        """
        Plots

        :param axs: DESCRIPTION, defaults to None
        :type axs: TYPE, optional
        :param fig: DESCRIPTION, defaults to None
        :type fig: TYPE, optional
        :param row: DESCRIPTION, defaults to 1
        :type row: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if fig is None:
            fig = plt.figure()

        plot_names = ['plot_axial_strain_x', 'plot_axial_strain_y', 'plot_shear_strain_xy']
        axs = []
        for i, name in enumerate(plot_names):
            axs.append(getattr(self, name)(ax=plt.subplot(row, 3, i + 1), fig=fig))

        if self.plane_stress:
            comment = 'Strain "Plane Stress"'
        else:
            comment = 'Strain "Plane Strain"'

        fig.suptitle(comment)

        return axs

    def plot_stress(self, axs=None, fig=None, row=1):
        """
        Plots

        :param axs: DESCRIPTION, defaults to None
        :type axs: TYPE, optional
        :param fig: DESCRIPTION, defaults to None
        :type fig: TYPE, optional
        :param row: DESCRIPTION, defaults to 1
        :type row: TYPE, optional

        :return: DESCRIPTION
        :rtype: TYPE
        """

        if fig is None:
            fig = plt.figure()

        plot_names = ['plot_axial_stress_x', 'plot_axial_stress_y', 'plot_shear_stress_xy']
        axs = []
        for i, name in enumerate(plot_names):
            axs.append(getattr(self, name)(ax=plt.subplot(row, 3, i + 1), fig=fig))

        if self.plane_stress:
            comment = 'Stress "Plane Stress"'
        else:
            comment = 'Stress "Plane Strain"'

        fig.suptitle(comment)

        return axs

    def shear_strain_xy(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        strain = self.strain
        shear_strain_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_strain_xy[element] = strain[element][2]

        return shear_strain_xy

    def shear_stress_xy(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        stress = self.stress
        shear_stress_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_stress_xy[element] = stress[element][2]

        return shear_stress_xy


class ElasticityResults3D(ElasticityResults):
    """
    This class
    """
    # _standalone_in_db = True
    # _non_serializable_attributes = []
    # _non_eq_attributes = ['name']
    # _non_hash_attributes = ['name']
    # _generic_eq = True
    # def __init__(self, mesh: vmmesh.Mesh, result_vector: List[float]):
    #     self.mesh = mesh
    #     self.result_vector = result_vector

    #     self.displacement_vectors_per_node = self._displacement_vectors_per_node()
    #     self.displacements_per_element = self._displacements_per_element()
    #     self.strain, self.stress = self._strain_stress_per_element()
    #     self.deformed_nodes = self._deformed_nodes()
    #     self.deformed_mesh = self._deformed_mesh()

    #     Result.__init__(self, mesh, result_vector)

    def axial_strain_x(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        strain = self.strain
        axial_strain_x = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_strain_x[element] = strain[element][0]

        return axial_strain_x

    def axial_strain_y(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        strain = self.strain
        axial_strain_y = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_strain_y[element] = strain[element][1]

        return axial_strain_y

    def axial_strain_z(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        strain = self.strain
        axial_strain_z = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_strain_z[element] = strain[element][2]

        return axial_strain_z

    def axial_stress_x(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        stress = self.stress
        axial_stress_x = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_stress_x[element] = stress[element][0]

        return axial_stress_x

    def axial_stress_y(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        stress = self.stress
        axial_stress_y = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_stress_y[element] = stress[element][1]

        return axial_stress_y

    def axial_stress_z(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        stress = self.stress
        axial_stress_z = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_stress_z[element] = stress[element][2]

        return axial_stress_z

    def displacement_per_node_z(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        return [displacement[2] for displacement in self.displacement_vectors_per_node]

    # def plot_axial_strain_z(self, ax=None, fig=None):

    #     return self.plot_constraints(constraint_name='axial_strain_z', ax=ax, fig=fig)

    # def plot_axial_stress_z(self, ax=None, fig=None):

    #     return self.plot_constraints(constraint_name='axial_stress_z', ax=ax, fig=fig)

    # def plot_displacement_per_node_z(self, ax=None):

    #     return self.plot_displacements(displacement_name='displacement_per_node_z', ax=ax)

    # def plot_shear_strain_yz(self, ax=None, fig=None):

    #     return self.plot_constraints(constraint_name='shear_strain_yz', ax=ax, fig=fig)

    # def plot_shear_strain_zx(self, ax=None, fig=None):

    #     return self.plot_constraints(constraint_name='shear_strain_zx', ax=ax, fig=fig)

    # def plot_shear_stress_yz(self, ax=None, fig=None):

    #     return self.plot_constraints(constraint_name='shear_stress_yz', ax=ax, fig=fig)

    # def plot_shear_stress_zx(self, ax=None, fig=None):

    #     return self.plot_constraints(constraint_name='shear_stress_zx', ax=ax, fig=fig)

    # def plot_strain(self, axs=None, fig=None, row=1):
    #     if fig is None:
    #         fig = plt.figure()

    #     plot_names = ['plot_axial_strain_x', 'plot_axial_strain_y', 'plot_shear_strain_xy',
    #                   'plot_axial_strain_z', 'plot_shear_strain_yz', 'plot_shear_strain_zx']
    #     axs = []
    #     for i, name in enumerate(plot_names):
    #         axs.append(getattr(self, name)(ax=plt.subplot(row, 3, i+1), fig=fig))

    #     return axs

    # def plot_stress(self, axs=None, fig=None, row=1):
    #     if fig is None:
    #         fig = plt.figure()

    #     plot_names = ['plot_axial_stress_x', 'plot_axial_stress_y', 'plot_shear_stress_xy',
    #                   'plot_axial_stress_z', 'plot_shear_stress_yz', 'plot_shear_stress_zx']
    #     axs = []
    #     for i, name in enumerate(plot_names):
    #         axs.append(getattr(self, name)(ax=plt.subplot(row, 3, i+1), fig=fig))

    #     return axs

    def shear_strain_xy(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        strain = self.strain
        shear_strain_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_strain_xy[element] = strain[element][3]

        return shear_strain_xy

    def shear_strain_yz(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        strain = self.strain
        shear_strain_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_strain_xy[element] = strain[element][4]

        return shear_strain_xy

    def shear_strain_zx(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        strain = self.strain
        shear_strain_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_strain_xy[element] = strain[element][5]

        return shear_strain_xy

    def shear_stress_xy(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        stress = self.stress
        shear_stress_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_stress_xy[element] = stress[element][3]

        return shear_stress_xy

    def shear_stress_yz(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        stress = self.stress
        shear_stress_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_stress_xy[element] = stress[element][4]

        return shear_stress_xy

    def shear_stress_zx(self):
        """
        Defines

        :return: DESCRIPTION
        :rtype: TYPE
        """

        stress = self.stress
        shear_stress_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_stress_xy[element] = stress[element][5]

        return shear_stress_xy
