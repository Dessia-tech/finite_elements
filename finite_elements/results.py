#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing objects related to finite elements analysis results
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as npy
import volmdlr as vm
import volmdlr.mesh as vmmesh
from dessia_common import DessiaObject
from typing import List  # Tuple, TypeVar
from finite_elements.core import MU, blue_red
import finite_elements.core
from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner
from matplotlib import cm


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

    @property
    def dimension(self):
        return self.mesh.elements_groups[0].elements[0].dimension


class MagneticResults(Result):
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

    def brbtetha(self, gap_elements_group):
        element_to_magnetic_field = self.magnetic_field_per_element
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

    def _magnetic_field_norm(self):
        element_to_magnetic_field = self.magnetic_field_per_element
        Bs = [B.norm() for B in list(element_to_magnetic_field.values())]

        return Bs

    def _magnetic_field_per_element(self):
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
                B_x = float(c1 *
                            self.result_vector[indexes[0]] +
                            c2 *
                            self.result_vector[indexes[1]] +
                            c3 *
                            self.result_vector[indexes[2]])
                B_y = float(-b1 * self.result_vector[indexes[0]] - b2 *
                            self.result_vector[indexes[1]] - b3 * self.result_vector[indexes[2]])
                element_to_magnetic_field[element] = vm.Vector2D(B_x, B_y)
        return element_to_magnetic_field

    def maxwell_stress_tensor(self, element):
        """
        Computes the Maxwell stress tensor for one element.
        Returns a list made of sigma_r_r, sigma_r_theta and sigma_theta_theta, \
        since sigma_r_theta = sigma_theta_r.

        :param element: The element on which the tensor is computed.
        :type element: a TriangularElement object
        """
        element_to_magnetic_field = self.magnetic_field_per_element
        vector_B = element_to_magnetic_field[element]
        element_center = element.center
        e_r = vm.Vector2D(element_center.vector)
        e_r.Normalize()
        e_teta = vm.Vector2D((-e_r[1], e_r[0]))
        B_r = vector_B.Dot(e_r)
        B_teta = vector_B.Dot(e_teta)

        sigma_rr = 1 / MU * B_r**2 - 1 / (2 * MU) * vector_B.Norm()**2
        sigma_rteta = 1 / MU * B_r * B_teta
        sigma_tetateta = 1 / MU * B_teta**2 - 1 / (2 * MU) * vector_B.Norm()**2
        sigma_rr_rteta_tetateta = [sigma_rr, sigma_rteta, sigma_tetateta]
        return sigma_rr_rteta_tetateta

    def torque(self, air_gap_elements_group_name, length_motor, radius_stator, radius_rotor):  # , nb_notches):
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
        element_to_magnetic_field = self.magnetic_field_per_element

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
        T = length_motor / (MU * (radius_stator - radius_rotor)) * somme

        return T

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
            color = jet(int(x * 256))[:3]

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

    def plot_magnetic_field(self, ax=None, Bmax=None):
        """
        Plots the mesh with colored triangular elements representing the \
        intensity of the induction field inside the Machine.
        """
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()
        element_to_magnetic_field = self.magnetic_field_per_element
        Bs = self.magnetic_field_norm

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

        element_to_magnetic_field = self.magnetic_field_per_element

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
        tri_refi, z_test_refi = refiner.refine_field(Z, subdiv=subdiv)

        levels = npy.arange(0., 1., 0.05)
        cmap = cm.get_cmap(name='Blues', lut=None)
        ax.tricontour(tri_refi, z_test_refi, levels=levels,  # cmap=cmap,
                      linewidths=[2.0, 0.5, 1.0, 0.5])
#        ax.triplot(tri_refi, color='0.97')
#        ax.triplot(tri, color='0.7')
#        ax.tricontour(x, y, Z)
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

        element_to_magnetic_field = self.magnetic_field_per_element
        Bs = self.magnetic_field_norm

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

    def plot_potential_vector(self, ax=None):
        """
        Plots the mesh with colored triangular elements representing the \
        intensity of the magnetic potential inside the Machine.
        """

        triang = finite_elements.core.get_triangulation(self.mesh)
        z = npy.asarray([p for p in self.result_vector[:len(self.mesh.nodes)]])  # p[0]
        z_min, z_max = min(z), max(z)

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()
        ax.tricontourf(triang, z, cmap=blue_red)
        ax.triplot(triang, 'k-')
        ax.set_title('Triangular grid')

        norm = mpl.colors.Normalize(vmin=z_min, vmax=z_max)
        sm = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ticks=npy.linspace(z_min, z_max, 10))
        cbar.set_label('Potential Vector')

        return ax


class ElasticityResults(Result):
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

        # self.displacement_vectors_per_node = self._displacement_vectors_per_node()
        # self.displacements_per_element = self._displacements_per_element()
        # self.energy_per_element = self._energy_per_element()
        # self.energy = self._energy()
        # self.strain, self.stress = self._strain_stress_per_element()
        # self.deformed_nodes = self._deformed_nodes()
        # self.deformed_mesh = self._deformed_mesh()

        self._displacement_vectors_per_node = None
        self._displacements_per_element = None
        self._energy_per_element = None
        self._energy = None
        self._strain, self._stress = None, None
        self._deformed_nodes = None
        self._deformed_mesh = None

        Result.__init__(self, mesh, result_vector)

    @property
    def deformed_mesh(self):
        if not self._deformed_mesh:
            self._deformed_mesh = self.deformed_mesh_m()
        return self._deformed_mesh

    def deformed_mesh_m(self, amplitude=1):
        if amplitude == 1:
            deformed_nodes = self.deformed_nodes
        else:
            deformed_nodes = self.deformed_nodes_m(amplitude=amplitude)

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

    @property
    def deformed_nodes(self):
        if not self._deformed_nodes:
            self._deformed_nodes = self.deformed_nodes_m()
        return self._deformed_nodes

    def deformed_nodes_m(self, amplitude=1):
        displacement_field_vectors = self.displacement_vectors_per_node
        deformed_nodes = []
        for i, node in enumerate(self.mesh.nodes):
            obj = getattr(vmmesh, f'Node{self.__class__.__name__[-2::]}')
            deformed_nodes.append(getattr(obj, 'from_point')(node + displacement_field_vectors[i] * amplitude))

        return deformed_nodes

    def displacement_per_node_x(self):

        return [displacement[0] for displacement in self.displacement_vectors_per_node]

    def displacement_per_node_y(self):

        return [displacement[1] for displacement in self.displacement_vectors_per_node]

    @property
    def displacement_vectors_per_node(self):
        if not self._displacement_vectors_per_node:
            self._displacement_vectors_per_node = self.displacement_vectors_per_node_m()
        return self._displacement_vectors_per_node

    def displacement_vectors_per_node_m(self):
        nodes_number = len(self.mesh.nodes)
        positions = finite_elements.core.global_matrix_positions(dimension=self.dimension,
                                                                 nodes_number=nodes_number)
        displacement_field_vectors = {}
        q = self.result_vector

        for n, node in enumerate(self.mesh.nodes):
            displacement = []
            for i in range(self.dimension):
                displacement.append(q[positions[(n, i + 1)]])

            displacement_field_vectors[node] = getattr(
                vm, f'Vector{self.__class__.__name__[-2::]}')(*displacement)
            # displacement_field_vectors.append(vm.Vector2D(*displacement))

        return displacement_field_vectors

    @property
    def displacements_per_element(self):
        if not self._displacements_per_element:
            self._displacements_per_element = self.displacements_per_element_m()
        return self._displacements_per_element

    def displacements_per_element_m(self):
        positions = finite_elements.core.global_matrix_positions(dimension=self.dimension,
                                                                 nodes_number=len(self.mesh.nodes))
        q = self.result_vector

        displacements_per_element = {}
        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:
                displacements = []
                indexes = [self.mesh.node_to_index[point] for point in element.points]
                for index in indexes:
                    for i in range(self.dimension):
                        displacements.append(q[positions[(index, i + 1)]])

                displacements_per_element[element] = displacements
                element.displacements = displacements

        return displacements_per_element

    @property
    def energy(self):
        if not self._energy:
            self._energy = self.energy_m()
        return self._energy

    def energy_m(self):
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

        return sum([value for value in self.energy_per_element.values()])

    @property
    def energy_per_element(self):
        if not self._energy_per_element:
            _ = self.displacements_per_element_m()
            self._energy_per_element = self.energy_per_element_m()
        return self._energy_per_element

    def energy_per_element_m(self):
        energy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                energy[element] = (element.energy(self.plane_strain, self.plane_stress))
        return energy

    @property
    def strain(self):
        if not self._strain:
            self._strain, self._stress = self.strain_stress_per_element_m()
        return self._strain

    @property
    def stress(self):
        if not self._stress:
            self._strain, self._stress = self.strain_stress_per_element_m()
        return self._stress

    def strain_stress_per_element_m(self):
        element_to_strain, element_to_stress = {}, {}
        for elements_group in self.mesh.elements_groups:
            for element in elements_group.elements:
                element_to_strain[element] = (npy.matmul(element.b_matrix, element.displacements))
                element.strain = element_to_strain[element]
                element_to_stress[element] = (npy.matmul(npy.matmul(element.d_matrix(plane_strain=self.plane_strain, plane_stress=self.plane_stress),
                                                                    element.b_matrix),
                                                         element.displacements))
                element.stress = element_to_stress[element]

        return element_to_strain, element_to_stress

    def generate_vtk_file(self, file_name):
        self.mesh._gmsh.to_vtk(file_name)
        with open(file_name) as f_in:
            with open(file_name+'_results', "w") as f_out:
                for line in f_in:
                    f_out.write(line)
        f_out.close()
        f_in.close()

        nodes_correction = self.mesh._nodes_correction
        displacements = []
        for node in self.mesh._gmsh.nodes[0]['all_nodes']:
            try:
                displacements.append(self.displacement_vectors_per_node[node])
            except KeyError:
                displacements.append(self.displacement_vectors_per_node[nodes_correction[node]])

        lines = ['POINT_DATA ' + str(len(self.mesh._gmsh.nodes[0]['all_nodes']))]
        lines.append('SCALARS ' + 'Displacement_Magnitude float 1')
        lines.append('LOOKUP_TABLE default')

        for displacement in displacements:
            lines.append(str(displacement.norm()))

        lines.append('VECTORS Displacement_Vectors float')
        if displacement.__class__.__name__[-2] == '2':
            for displacement in displacements:
                lines.append(str([*displacement])[1:-1].replace(',','')+ ' 0')
        else:
            for displacement in displacements:
                lines.append(str([*displacement])[1:-1].replace(',',''))

        with open(file_name+'_results', "a+") as f_out:
            for line in lines:
                f_out.write(line)
                f_out.write('\n')
        f_out.close()

    '''
    def update_vtk_with_results(self, input_file_name, output_file_name):
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
        '''

class ElasticityResults2D(ElasticityResults):
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

        strain = self.strain
        axial_strain_x = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_strain_x[element] = strain[element][0]

        return axial_strain_x

    def axial_strain_y(self):

        strain = self.strain
        axial_strain_y = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_strain_y[element] = strain[element][1]

        return axial_strain_y

    def axial_stress_x(self):

        stress = self.stress
        axial_stress_x = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_stress_x[element] = stress[element][0]

        return axial_stress_x

    def axial_stress_y(self):

        stress = self.stress
        axial_stress_y = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_stress_y[element] = stress[element][1]

        return axial_stress_y

    def plot_axial_strain_x(self, ax=None, fig=None):

        return self.plot_constraints(constraint_name='axial_strain_x', ax=ax, fig=fig)

    def plot_axial_strain_y(self, ax=None, fig=None):

        return self.plot_constraints(constraint_name='axial_strain_y', ax=ax, fig=fig)

    def plot_axial_stress_x(self, ax=None, fig=None):

        return self.plot_constraints(constraint_name='axial_stress_x', ax=ax, fig=fig)

    def plot_axial_stress_y(self, ax=None, fig=None):

        return self.plot_constraints(constraint_name='axial_stress_y', ax=ax, fig=fig)

    def plot_constraints(self, constraint_name: str, ax=None, fig=None):

        if ax is None:
            fig, ax = plt.subplots()

        if hasattr(self, constraint_name):
            result = getattr(self, constraint_name)()
            result_values = [value for value in result.values()]
        else:
            raise NotImplementedError(
                f'Class {self.__class__.__name__} does not implement {constraint_name}')

        deformed_mesh = self.deformed_mesh
        B_max, B_min = finite_elements.core.get_bmin_bmax(result_values, Bmin=None, Bmax=None)
        B_to_color = finite_elements.core.get_colors(result_values, B_max=B_max, B_min=B_min)
        for g, group in enumerate(deformed_mesh.elements_groups):
            for e, element in enumerate(group.elements):
                element.plot(ax=ax, color=B_to_color[result[self.mesh.elements_groups[g].elements[e]]], fill=True)

        norm = mpl.colors.Normalize(vmin=B_min, vmax=B_max)
        sm = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ticks=npy.linspace(B_min, B_max, 10), ax=ax)
        cbar.set_label(constraint_name)

        ax.set_xlabel('x')
        ax.set_ylabel('y')

        return ax

    def plot_deformed_mesh(self, ax=None, amplitude=1):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        if amplitude != 1:
            self.deformed_mesh_m(amplitude=amplitude).plot(ax=ax)
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

        if hasattr(self, displacement_name):
            x = getattr(self, displacement_name)()
        elif displacement_name == 'displacement_per_node_xy':
            x = [displacement.norm() for displacement in self.displacement_vectors_per_node]
            displacement_name = 'displacement_magnitude_norm'
        else:
            raise NotImplementedError(
                f'Class {self.__class__.__name__} does not implement {displacement_name}')

        if amplitude != 1:
            mesh_fe = self.deformed_mesh_m(amplitude=amplitude)
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
        sm = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ticks=npy.linspace(x_min, x_max, 10))
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

        return self.plot_displacements(displacement_name='displacement_per_node_x', ax=ax, amplitude=amplitude)

    def plot_displacement_per_node_xy(self, ax=None, amplitude=1):

        return self.plot_displacements(displacement_name='displacement_per_node_xy', ax=ax, amplitude=amplitude)

    def plot_displacement_per_node_y(self, ax=None, amplitude=1):

        return self.plot_displacements(displacement_name='displacement_per_node_y', ax=ax, amplitude=amplitude)

    def plot_displacement_vectors_per_node(self, ax=None, amplitude=0.05):
        if ax is None:
            fig, ax = plt.subplots()

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

        result = self.energy_per_element
        result_values = [value for value in self.energy_per_element.values()]

        if amplitude != 1:
            deformed_mesh = self.deformed_mesh_m(amplitude=amplitude)
        else:
            deformed_mesh = self.deformed_mesh

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()

        B_max, B_min = finite_elements.core.get_bmin_bmax(result_values, Bmin=None, Bmax=None)
        B_to_color = finite_elements.core.get_colors(result_values, B_max=B_max, B_min=B_min)
        for g, group in enumerate(deformed_mesh.elements_groups):
            for e, element in enumerate(group.elements):
                element.plot(ax=ax, color=B_to_color[result[self.mesh.elements_groups[g].elements[e]]], fill=True)

        norm = mpl.colors.Normalize(vmin=B_min, vmax=B_max)
        sm = plt.cm.ScalarMappable(cmap=blue_red, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ticks=npy.linspace(B_min, B_max, 10), ax=ax)
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

        return self.plot_constraints(constraint_name='shear_strain_xy', ax=ax, fig=fig)

    def plot_shear_stress_xy(self, ax=None, fig=None):

        return self.plot_constraints(constraint_name='shear_stress_xy', ax=ax, fig=fig)

    def plot_strain(self, axs=None, fig=None, row=1):
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

        strain = self.strain
        shear_strain_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_strain_xy[element] = strain[element][2]

        return shear_strain_xy

    def shear_stress_xy(self):

        stress = self.stress
        shear_stress_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_stress_xy[element] = stress[element][2]

        return shear_stress_xy


class ElasticityResults3D(ElasticityResults):
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

        strain = self.strain
        axial_strain_x = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_strain_x[element] = strain[element][0]

        return axial_strain_x

    def axial_strain_y(self):

        strain = self.strain
        axial_strain_y = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_strain_y[element] = strain[element][1]

        return axial_strain_y

    def axial_strain_z(self):

        strain = self.strain
        axial_strain_z = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_strain_z[element] = strain[element][2]

        return axial_strain_z

    def axial_stress_x(self):

        stress = self.stress
        axial_stress_x = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_stress_x[element] = stress[element][0]

        return axial_stress_x

    def axial_stress_y(self):

        stress = self.stress
        axial_stress_y = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_stress_y[element] = stress[element][1]

        return axial_stress_y

    def axial_stress_z(self):

        stress = self.stress
        axial_stress_z = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                axial_stress_z[element] = stress[element][2]

        return axial_stress_z

    def displacement_per_node_z(self):

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

        strain = self.strain
        shear_strain_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_strain_xy[element] = strain[element][3]

        return shear_strain_xy

    def shear_strain_yz(self):

        strain = self.strain
        shear_strain_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_strain_xy[element] = strain[element][4]

        return shear_strain_xy

    def shear_strain_zx(self):

        strain = self.strain
        shear_strain_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_strain_xy[element] = strain[element][5]

        return shear_strain_xy

    def shear_stress_xy(self):

        stress = self.stress
        shear_stress_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_stress_xy[element] = stress[element][3]

        return shear_stress_xy

    def shear_stress_yz(self):

        stress = self.stress
        shear_stress_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_stress_xy[element] = stress[element][4]

        return shear_stress_xy

    def shear_stress_zx(self):

        stress = self.stress
        shear_stress_xy = {}
        for group in self.mesh.elements_groups:
            for element in group.elements:
                shear_stress_xy[element] = stress[element][5]

        return shear_stress_xy
