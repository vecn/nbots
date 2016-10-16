# Numerical Bots

## Abstract

Numerical bots is an incubating project which aims to implement the most used numerical calculations in scientific computing focusing on achieving the best performance (processor and memory efficiency) for single and parallel arquitectures, such as shared memory, distributed memory and general purpose graphic processors.

## Introduction

The project nbots is a self-contained ecosystem of numerical bots.
A numerical bot is an interdependent module which perform several related numerical tasks.

### List of bots
- Memory bot
- Container bot
- Math bot
- Graph bot
- Geometric bot
- Graphics bot
- Image bot
- Solver bot
- Optimization bot*
- Metaheuristics bot*
- Statistics bot
- PDE bot
- Topopt bot*

The '*' is for those bots which are not implemented yet.

## Compilation

Run 'gradle assemble' or use the shell script './gradlew assemble' (gradlew.bat for Windows) if gradle is not installed in your system.

## Acknowledges

The Center of Research in Mathematics, CIMAT, in many ways has been the backbone of this project, especially **the people** forming the department of Computational Sciences and its post-graduate programs, a big **thank you** to all of them.

## Developers and contributions

- Victor E. Cardoso (victorc@cimat.mx)
     * Library arch.
     * Bots integration and Unit tests.
     * Containers.
     * 2D Geometric algorithms (meshes, alpha-shapes, Voronoi, etc).
     * Finite element assemblers.
     * Linear algebra routines,
     * AMD: Graph labeling routines.
     * Spectral bissection,
     * Graphics bot (rasterizers).
- Miguel Vargas Felix (miguelvargas@cimat.mx)
     * System macros.
     * Parallel sparse solvers and structures.
     * Domain segmentation routines.
     * Nested dissections (graph labeling routine)
- Jorge Lopez (jorge.lopez@cimat.mx)
     * 3D tetrahedral mesher based on octree.
- Dora Elisa Alvarado (dora.alvarado@cimat.mx)
     * Image processing routines and filters.
- Miguel Angel Ochoa (miguel.ochoa@cimat.mx)
     * SIMP Topology optimization method.
- J. Gerardo Fuentes
     * Finite element electromagnetic solver.
- Sean T. Barrett (https://github.com/nothings/stb)
     * Pixmap exporters to JPG and PNG.
    