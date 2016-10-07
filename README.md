Numerical Bots
==============

Abstract
--------
Numerical bots is an incubating project which aims to implement the most used numerical calculations in scientific computing focusing on achieving the best performance (processor and memory efficiency) for single and parallel arquitectures, such as shared memory, distributed memory and general purpose graphic processors.

Introduction
------------
A numerical bot is an automatic tool to perform several related numerical tasks.
For example, the Geometric Bot calculates Delaunay triangulations, Finite Element Meshes, Voronoi diagrams, Alpha-shapes, Convex-Hulls, etc. and the Eigen Bot is a linear algebra library, which includes special routines for sparse systems.

The project nbots is a self-contained ecosystem of numerical bots.

Compilation
-----------
Run 'gradle assemble' or use the shell script './gradlew assemble' (gradlew.bat for Windows) if gradle is not installed in your system.

Developers and contributions
-----------------------
- Victor E. Cardoso (victorc@cimat.mx)
     > Library arch, 2D meshes, alpha-shapes, containers,
       finite element assemblers, linear algebra routines,
       graph labeling routines, spectral bissection,
       and graphics functions (rasterizers).
     > Bots integration and Unit tests.
- Miguel Vargas Felix (miguelvargas@cimat.mx)
     > Parallel sparse solvers and structures.
     > Domain segmentation routines.
- Jorge Lopez (jorge.lopez@cimat.mx)
     > 3D tetrahedral mesher based on octree.
- Dora Elisa Alvarado (dora.alvarado@cimat.mx)
     > Image processing routines and filters.
- Miguel Angel Ochoa (miguel.ochoa@cimat.mx)
     > SIMP Topology optimization method.
- Sean T. Barrett (https://github.com/nothings/stb)
     > Pixmap exporters to JPG and PNG.
    