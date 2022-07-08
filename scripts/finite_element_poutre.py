#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 16:18:29 2020
@author: ringhausen
"""
import finite_elements as fe
import finite_elements.elements
import finite_elements.loads
import finite_elements.analysis
import volmdlr.mesh as mesh 
import volmdlr as vm
import math
phase1 = 3
phase2 = 3
phase3 = 3
mu1 = 4*math.pi*1e-7 * 100000 # IRON
mu2 = 4*math.pi*1e-7          # VOID
mu3 = 4*math.pi*1e-7 * 50000  # IRON

elements_phase1 = []
for i in range(phase1):
    element1 = mesh.TriangularElement2D([vm.Point2D(i,0), vm.Point2D(i+1,0), vm.Point2D(i,1)])
    element2 = mesh.TriangularElement2D([vm.Point2D(i+1,1), vm.Point2D(i+1,0), vm.Point2D(i,1)])
    elements_phase1.extend([element1, element2])

magnetic_elements = [fe.elements.MagneticElement2D(element, mu1) for element in elements_phase1]
# group_phase1 = fe.elements.MagneticElementsGroup(magnetic_elements, mu1, 'phase1')
group_phase1 = mesh.ElementsGroup(magnetic_elements, 'phase1')

elements_phase2 = []
for i in range(phase1, phase1+phase2):
    element1 = mesh.TriangularElement2D([vm.Point2D(i,0), vm.Point2D(i+1,0), vm.Point2D(i,1)])
    element2 = mesh.TriangularElement2D([vm.Point2D(i+1,1), vm.Point2D(i+1,0), vm.Point2D(i,1)])
    elements_phase2.extend([element1, element2])

magnetic_elements = [fe.elements.MagneticElement2D(element, mu2) for element in elements_phase2]
# group_phase2 = fe.elements.MagneticElementsGroup(magnetic_elements, mu2, 'phase2')
group_phase2 = mesh.ElementsGroup(magnetic_elements, 'phase2')


elements_phase3 = []
for i in range(phase1+phase2, phase1+phase2+phase3):
    element1 = mesh.TriangularElement2D([vm.Point2D(i,0), vm.Point2D(i+1,0), vm.Point2D(i,1)])
    element2 = mesh.TriangularElement2D([vm.Point2D(i+1,1), vm.Point2D(i+1,0), vm.Point2D(i,1)])
    elements_phase3.extend([element1, element2])

magnetic_elements = [fe.elements.MagneticElement2D(element, mu3) for element in elements_phase3]
# group_phase3 = fe.elements.MagneticElementsGroup(magnetic_elements, mu3, 'phase3')
group_phase3 = mesh.ElementsGroup(magnetic_elements, 'phase3')

mesh = mesh.Mesh([group_phase1, group_phase2, group_phase3])
elements_loads = [fe.loads.ElementsLoad(group_phase1.elements[0:2], 1e10, 1)]
node_loads = []
for node in mesh.nodes:
    if math.isclose(node[0], phase1+phase2+phase3, abs_tol=1e-6):
        node_loads.append(fe.loads.NodeLoad(node, 0, 1))
analysis = fe.analysis.FiniteElementAnalysis(mesh, elements_loads, node_loads, [], [], [], [])
analysis.plot_elements_loads()
analysis.plot_elements_permeability()
results = analysis.solve()
results.plot_magnetic_field()