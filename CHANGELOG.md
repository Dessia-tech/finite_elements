# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.1.0] - unreleased

### Refactorings

* Divide electromag.py into different small thematic files
  elements.py, conditions.py, loads.py, analysis.py, results.py, core.py
* Generalize methods to be used in different fields (magnetic, elasticity, etc.)
* Update classes and methods:
   - elementary_source_matrix > element_to_node_factors
   - ConstantLoad > ElementsLoad
   - SingleNodeLoad > NodeLoad

### New Features

* Add new classes and methods:
   - MagneticElement2D
   - BoundaryCondition, NodeBoundaryCondition, EdgeBoundaryCondition, ElementBoundaryCondition
   - ElementLoad, EdgeLoad
   - FiniteElementAnalysis: modal_analysis
   - MagneticResults
   - ElasticityElement: d_matrix, energy
   - ElasticityTriangularElement2D: d_matrix(plane_strain, plane_stress), b_matrix, k_matrix, strain, stress, m_matrix
   - ElasticityTetrahedralElement3D: d_matrix, b_matrix, k_matrix, m_matrix
   - ElasticityResults: displacements, stress, strain, deformed mesh, energy, vtk file
   - ElasticityResults2D: axial_strain/stress, shear_strain/stress, different plots
   - ElasticityResults3D: axial_strain/stress, shear_strain/stress
   - Material
* Add new scripts with usecases and examples for tests (Magnetic, Elasticity2D)

### Performance improvements

* Update ElasticityResults (use property and hidden attributes)
*

### Fixed

*
*


## [v0.0.1] - 31/08/2020
