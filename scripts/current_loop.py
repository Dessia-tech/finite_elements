#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 17:39:02 2020
@author: ringhausen
"""
# import electric_machines.core as em
import volmdlr as vm
import volmdlr.mesh as vmmesh
import matplotlib.pyplot as plt 
import math
import finite_elements
import finite_elements as fe
import finite_elements.elements
import finite_elements.loads
import finite_elements.elements_analysis

MU = 4*math.pi*1e-7

# %% Method

def create_loop(vertical_resolution, horizontal_resolution, mu_r_loop, box_multiplier, offset_radius):
    
############################## LEFT WIRE ##############################
    
    pt_low_left = vm.Point2D(-0.375, -0.125/2) + vm.Point2D(-0.125*offset_radius, 0)
    pt_low_right = vm.Point2D(-0.25, -0.125/2) + vm.Point2D(-0.125*offset_radius, 0)
    pt_up_left = vm.Point2D(-0.375, 0.125/2) + vm.Point2D(-0.125*offset_radius, 0)
#    pt_up_right = vm.Point2D(-0.25, 0.125)
    
    low_points = [pt_low_left]
    low_segment = vm.edges.LineSegment2D(pt_low_left, pt_low_right)
    for i in range(1, horizontal_resolution):
        new_low_pt = low_segment.point_at_abscissa(i * low_segment.length() / horizontal_resolution)
        low_points.append(new_low_pt)
    low_points.append(pt_low_right)
    
    left_points = [pt_low_left]
    left_segment = vm.edges.LineSegment2D(pt_low_left, pt_up_left)
    for i in range(1, vertical_resolution):
        new_left_pt = left_segment.point_at_abscissa(i * left_segment.length() / vertical_resolution)
        left_points.append(new_left_pt)
    left_points.append(pt_up_left)
    
    all_points_left = []
    for pt_left in left_points:
        for pt_low in low_points:
            new_pt = vm.Point2D(pt_low[0], pt_left[1])
            all_points_left.append(new_pt)
            
    elements = []
    for i in range(len(all_points_left)-horizontal_resolution-2):
        if (i+1) % (horizontal_resolution+1) == 0 and i != 0:
            continue
        elem1 = vmmesh.TriangularElement2D([all_points_left[i], all_points_left[i+1], all_points_left[i+horizontal_resolution+2]])
        elem2 = vmmesh.TriangularElement2D([all_points_left[i], all_points_left[i+horizontal_resolution+1], all_points_left[i+horizontal_resolution+2]])
        elements.extend([elem1, elem2])
        
    left_elem_group = fe.elements.MagneticElementsGroup(elements, mu_r_loop*MU, 'Left wire')
############################## RIGHT WIRE ##############################    
    
    pt_low_left = vm.Point2D(0.25, -0.125/2) + vm.Point2D(0.125*offset_radius, 0)
    pt_low_right = vm.Point2D(0.375, -0.125/2) + vm.Point2D(0.125*offset_radius, 0)
    pt_up_left = vm.Point2D(0.25, 0.125/2) + vm.Point2D(0.125*offset_radius, 0)
#    pt_up_right = vm.Point2D(-0.25, 0.125)
    
    low_points = [pt_low_left]
    low_segment = vm.edges.LineSegment2D(pt_low_left, pt_low_right)
    for i in range(1, horizontal_resolution):
        new_low_pt = low_segment.point_at_abscissa(i * low_segment.length() / horizontal_resolution)
        low_points.append(new_low_pt)
    low_points.append(pt_low_right)
    
    left_points = [pt_low_left]
    left_segment = vm.edges.LineSegment2D(pt_low_left, pt_up_left)
    for i in range(1, vertical_resolution):
        new_left_pt = left_segment.point_at_abscissa(i * left_segment.length() / vertical_resolution)
        left_points.append(new_left_pt)
    left_points.append(pt_up_left)
    
    all_points_right = []
    for pt_left in left_points:
        for pt_low in low_points:
            new_pt = vm.Point2D(pt_low[0], pt_left[1])
            all_points_right.append(new_pt)
            
    elements = []
    for i in range(len(all_points_right)-horizontal_resolution-2):
        if (i+1) % (horizontal_resolution+1) == 0 and i != 0:
            continue
        elem1 = vmmesh.TriangularElement2D([all_points_right[i], all_points_right[i+1], all_points_right[i+horizontal_resolution+2]])
        elem2 = vmmesh.TriangularElement2D([all_points_right[i], all_points_right[i+horizontal_resolution+1], all_points_right[i+horizontal_resolution+2]])
        elements.extend([elem1, elem2])
        
    right_elem_group = fe.elements.MagneticElementsGroup(elements, mu_r_loop*MU, 'Right wire')
    
############################## BOX ##############################
    
    horizontal_resolution = horizontal_resolution * 4 * 4   
    vertical_resolution = vertical_resolution * 8 * 2
    
    pt_low_left = vm.Point2D(-1, -1) * box_multiplier
    pt_low_right = vm.Point2D(1, -1) * box_multiplier
    pt_up_left = vm.Point2D(-1, 1) * box_multiplier
    
    low_points = [pt_low_left]
    low_segment = vm.edges.LineSegment2D(pt_low_left, pt_low_right)
    for i in range(1, horizontal_resolution):
        new_low_pt = low_segment.point_at_abscissa(i * low_segment.length() / horizontal_resolution)
        low_points.append(new_low_pt)
    low_points.append(pt_low_right)
    
    left_points = [pt_low_left]
    left_segment = vm.edges.LineSegment2D(pt_low_left, pt_up_left)
    for i in range(1, vertical_resolution):
        new_left_pt = left_segment.point_at_abscissa(i * left_segment.length() / vertical_resolution)
        left_points.append(new_left_pt)
    left_points.append(pt_up_left)
    
    all_points = []
    for pt_left in left_points:
        for pt_low in low_points:
            new_pt = vm.Point2D(pt_low[0], pt_left[1])
            all_points.append(new_pt)
            
    elements = []
    for i in range(len(all_points)-horizontal_resolution-2):
        if (i+1) % (horizontal_resolution+1) == 0 and i != 0:
            continue
        elif all_points[i] in all_points_left[:-int(horizontal_resolution/16)-1] \
        and not math.isclose(all_points[i][0], -0.25 - 0.125*offset_radius, abs_tol=1e-6):
            continue
        elif all_points[i] in all_points_right[:-int(horizontal_resolution/16)-1] \
        and not math.isclose(all_points[i][0], 0.25+0.125 + 0.125*offset_radius, abs_tol=1e-6):
            continue
        elem1 = vmmesh.TriangularElement2D([all_points[i], all_points[i+1], all_points[i+horizontal_resolution+2]])
        elem2 = vmmesh.TriangularElement2D([all_points[i], all_points[i+horizontal_resolution+1], all_points[i+horizontal_resolution+2]])
        elements.extend([elem1, elem2])
        
    box_elem_group = fe.elements.MagneticElementsGroup(elements, MU, 'Box')
    
    return vmmesh.Mesh([left_elem_group, right_elem_group, box_elem_group])

# %% Example

# Set intensity
intensity = 1000
box_multiplier = 3
offset_radius = 5
# The resolutions must be even numbers
#mesh = create_loop(4, 4, em.Steel.relative_permeability)
mesh = create_loop(4, 4, 1, box_multiplier, offset_radius)
    
node_loads = []
for node in mesh.nodes:
    if math.isclose(node[0], -1 * box_multiplier, abs_tol=1e-6) \
    or math.isclose(node[0], 1 * box_multiplier, abs_tol=1e-6) \
    or math.isclose(node[1], -1 * box_multiplier, abs_tol=1e-6) \
    or math.isclose(node[1], 1 * box_multiplier, abs_tol=1e-6):
        node_loads.append(fe.loads.SingleNodeLoad(node, 0))
element_load = [fe.loads.ConstantLoad(mesh.elements_groups[0].elements, intensity),
                fe.loads.ConstantLoad(mesh.elements_groups[1].elements, -intensity)]
analysis = fe.elements_analysis.FiniteElementAnalysis(mesh, element_load, node_loads, [], [])
result = analysis.solve()
ax = mesh.elements_groups[0].plot()
mesh.elements_groups[1].plot(ax=ax)
result.plot_magnetic_field_vectors(ax=ax, amplitude=0.05) #, Bmax=0.00003)