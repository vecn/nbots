# Developers and contributions

## Numerical bots

|ID |          Name           |                Caption                 |
|---|-------------------------|----------------------------------------|
|mem| Memory bot              | Handles memory allocations             |
|cnt| Container bot           | Process dynamic collections            |
|mat| Math bot                | Standard numerical procedures          |
|g2d| Graphics bot            | 2D-vector drawing functions            |
|grp| Graph bot               | Graph-specific routines                |
|sol| Solver bot              | Linear algebra operations              |
|opt| Optimization bot        | Constrained and unconstrained optim    |
|sta| Statistics bot          | Pseudo-random number generators        |
|geo| Geometric bot           | Implements geometric algorithms        |
|img| Image bot               | Image processing routines and filters  |
|mhe| Metaheuristics bot      | Black-box and combinatorial optim      |
|pde| PDE bot                 | Solve PDE using FEM & CVFA             |
|top| Topopt bot              | Topology optimization problems         |


## Developers

| ID |       Name        |      github        |
|----|-------------------|--------------------|
|vic |Victor E. Cardoso  | @vecn              |
|mvf |Miguel Vargas-Felix| @MiguelVargasFelix |
|jol |Jorge Lopez        | @mcjorgelopez      |
|dal |Dora E. Alvarado   | @alvaradocde       |
|mik |Miguel Angel Ochoa | @miguel-ochoa      |
|ger |J. Gerardo Fuentes | @Almeida88         |
|eot |Ernesto Ortega     | @netohomes         |

## Contributions

| Dev / Bot | mem | cnt | mat | g2d | grp | sol | opt | sta | geo | img | mhe | pde | top |
|:---------:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|  **vic**  |  x  |  x  |  x  |  x  |  x  |  x  |     |  x  |  x  |  x  |     |  x  |     |
|  **mvf**  |  x  |     |  x  |     |  x  |  x  |     |     |     |     |     |  x  |     |
|  **jol**  |     |     |     |     |     |     |     |     |  x  |     |     |     |     |
|  **dal**  |     |     |     |     |     |     |     |     |     |  x  |     |     |     |
|  **mik**  |     |     |     |     |     |     |  x  |     |     |     |     |     |  x  |
|  **ger**  |     |     |     |     |     |     |     |     |     |     |     |  x  |     |

## Detailed contributions

| Project task                                | Bot |  Developers  |
|---------------------------------------------|-----|--------------|
|Bots integration                             | all |vic           |
|System macros                                | all |mvf           |
|Documentation                                | all |vic, mvf, eot |
|Memory allocations                           | mem |vic, mvf      |
|Membank for structs of same length           | mem |vic           |
|Hash table                                   | cnt |vic           |
|Queue and Stack (Linked list)                | cnt |vic           |
|Sorted struct (AVL Tree)                     | cnt |vic           |
|Heap (Half sorted tree)                      | cnt |vic           |
|Array procedures (qsort, bsearch, swap, etc) | cnt |vic           |
|Generic Iterators                            | cnt |vic           |      
|2D Delaunay                                  | geo |vic           |
|2D Alpha-complex                             | geo |vic           |
|2D meshes                                    | geo |vic           |
|Static elasticity FEM assembler              | pde |vic, mvf      |
|Approximated minimum degree                  | grp |vic           |
|Spectral bisection                           | grp |vic           |
|Rasterizer                                   | g2d |vic           |
|Par parse solvers and structures             | sol |mvf, vic      |
|Domain segmentation                          | grp |mvf           |
|Nested dissection                            | grp |mvf           |
|3D Tetrahedral mesher (octree)               | geo |jol           |

## Contributions from 'Public domain code'
- (img) Image reader and writers (PNG, JPG, BMP).
  Sean T. Barrett
  ([GitHub repo](https://github.com/nothings/stb))
- (geo) Robust predicates.
  Jonathan Richard Shewchuk
  ([Web page](https://www.cs.cmu.edu/~quake/robust.html))
- (grp) Nested Dissection.
  George Karypis.
  ([Web page](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview))
    