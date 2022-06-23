#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 18:08:16 2020
@author: ringhausen
"""
import electric_machines.core as em
import volmdlr as vm
import volmdlr.mesh as vmmesh
import math
import electric_machines.electromag_finite_elements as fe
MU = 4*math.pi*1e-7
def create_magnet(horizontal_resolution, vertical_resolution, mu_r_magnet):
    
    pt_low_left = vm.Point2D((-0.25, -0.125))
    pt_low_right = vm.Point2D((0.25, -0.125))
    pt_up_left = vm.Point2D((-0.25, 0.125))
    
    low_points = [pt_low_left]
    low_segment = vm.LineSegment2D(pt_low_left, pt_low_right)
    for i in range(1, horizontal_resolution):
        new_low_pt = low_segment.PointAtCurvilinearAbscissa(i * low_segment.Length() / horizontal_resolution)
        low_points.append(new_low_pt)
    low_points.append(pt_low_right)
    
    left_points = [pt_low_left]
    left_segment = vm.LineSegment2D(pt_low_left, pt_up_left)
    for i in range(1, vertical_resolution):
        new_left_pt = left_segment.PointAtCurvilinearAbscissa(i * left_segment.Length() / vertical_resolution)
        left_points.append(new_left_pt)
    left_points.append(pt_up_left)
    
    all_points_manget = []
    for pt_left in left_points:
        for pt_low in low_points:
            new_pt = vm.Point2D((pt_low[0], pt_left[1]))
            all_points_manget.append(new_pt)
            
    elements = []
    for i in range(len(all_points_manget)-horizontal_resolution-2):
        if (i+1) % (horizontal_resolution+1) == 0 and i != 0:
            continue
        elem1 = vmmesh.TriangularElement([all_points_manget[i], all_points_manget[i+1], all_points_manget[i+horizontal_resolution+2]])
        elem2 = vmmesh.TriangularElement([all_points_manget[i], all_points_manget[i+horizontal_resolution+1], all_points_manget[i+horizontal_resolution+2]])
        elements.extend([elem1, elem2])
        
    magnet_elem_group = vmmesh.ElementsGroup(elements, mu_r_magnet*MU, 'Magnet')
    
##############################################################################
    
    horizontal_resolution = horizontal_resolution * 4
    vertical_resolution = vertical_resolution * 8
    
    pt_low_left = vm.Point2D((-1, -1))
    pt_low_right = vm.Point2D((1, -1))
    pt_up_left = vm.Point2D((-1, 1))
    
    low_points = [pt_low_left]
    low_segment = vm.LineSegment2D(pt_low_left, pt_low_right)
    for i in range(1, horizontal_resolution):
        new_low_pt = low_segment.PointAtCurvilinearAbscissa(i * low_segment.Length() / horizontal_resolution)
        low_points.append(new_low_pt)
    low_points.append(pt_low_right)
    
    left_points = [pt_low_left]
    left_segment = vm.LineSegment2D(pt_low_left, pt_up_left)
    for i in range(1, vertical_resolution):
        new_left_pt = left_segment.PointAtCurvilinearAbscissa(i * left_segment.Length() / vertical_resolution)
        left_points.append(new_left_pt)
    left_points.append(pt_up_left)
    
    all_points = []
    for pt_left in left_points:
        for pt_low in low_points:
            new_pt = vm.Point2D((pt_low[0], pt_left[1]))
            all_points.append(new_pt)
            
    elements = []
    for i in range(len(all_points)-horizontal_resolution-2):
        if (i+1) % (horizontal_resolution+1) == 0 and i != 0:
            continue
        elif all_points[i] in all_points_manget[:-int(horizontal_resolution/4)-1] \
        and not math.isclose(all_points[i][0], 0.25, abs_tol=1e-6):
            continue
        elem1 = vmmesh.TriangularElement([all_points[i], all_points[i+1], all_points[i+horizontal_resolution+2]])
        elem2 = vmmesh.TriangularElement([all_points[i], all_points[i+horizontal_resolution+1], all_points[i+horizontal_resolution+2]])
        elements.extend([elem1, elem2])
        
    box_elem_group = vmmesh.ElementsGroup(elements, MU, 'Box')
    
    return vmmesh.Mesh([magnet_elem_group, box_elem_group])
# magnet = em.RECOMA26
magnet = Magnet(1.04, 765000, 8300)
# The resolutions must be even numbers
mesh = create_magnet(6, 6, magnet.relative_permeability)
node_loads = []
for node in mesh.nodes:
    if math.isclose(node[0], -1, abs_tol=1e-6) \
    or math.isclose(node[0], 1, abs_tol=1e-6) \
    or math.isclose(node[1], -1, abs_tol=1e-6) \
    or math.isclose(node[1], 1, abs_tol=1e-6):
        node_loads.append(fe.SingleNodeLoad(node, 0))
magnet_loads = [fe.MagnetLoad(mesh.elements_groups[0].elements, [], -1/MU * magnet.remanent_induction * vm.Y2D)]
analysis = fe.FiniteElementAnalysis(mesh, [], node_loads, magnet_loads, [])
result = analysis.solve()
ax = mesh.elements_groups[0].plot()
result.plot_magnetic_field_vectors(ax=ax, amplitude=0.05)
ax = mesh.elements_groups[0].plot()
result.plot_magnetic_field_contour(ax=ax)
