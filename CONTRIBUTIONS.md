# Developers and contributions

## Numerical bots

|ID |          Name           |                Caption                 |
|---|-------------------------|----------------------------------------|
|mem| Memory bot              | Handles memory allocations             |
|cnt| Container bot           | Process dynamic collections            |
|mat| Math bot                | Standard numerical procedures          |
|g2d| Graphics bot            | 2D-vector drawing functions            |
|grp| Graph bot               | Graph-specific routines                |
|geo| Geometric bot           | Implements geometric algorithms        |
|img| Image bot               | Image processing routines and filters  |
|sol| Solver bot              | Linear algebra operations              |
|opt| Optimization bot        | Constrained and unconstrained optim    |
|mhe| Metaheuristics bot      | Black-box and combinatorial optim      |
|sta| Statistics bot          | Pseudo-random number generators        |
|pde| PDE bot                 | Solve PDE using FEM & CVFA             |
|top| Topopt bot              | Topology optimization problems         |


## Developers

| ID |       Name        |      github        |          email        |
|----|-------------------|--------------------|-----------------------|
|vic |Victor E. Cardoso  | @vecn              | victorc@cimat.mx      |
|mvf |Miguel Vargas-Felix| @MiguelVargasFelix | miguelvargas@cimat.mx |
|jol |Jorge Lopez        | @j1                | jorge.lopez@cimat.mx  |
|dal |Dora E. Alvarado   | @d1                | dora.alvarado@cimat.mx|
|mik |Miguel Angel Ochoa | @mo                | miguel.ochoa@cimat.mx |
|ger |J. Gerardo Fuentes | @jg                | juan.fuentes@cimat.mx |

## Contributions

| Dev / Bot | mem | cnt | mat | g2d | grp | geo | img | sol | opt | mhe | sta | pde | top |
|:---------:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|  **vic**  |  x  |  x  |  x  |  x  |  x  |  x  |  x  |  x  |     |     |  x  |  x  |     |
|  **mvf**  |  x  |     |  x  |     |  x  |     |     |  x  |     |     |     |  x  |     |
|  **jol**  |     |     |     |     |     |  x  |     |     |     |     |     |     |     |
|  **dal**  |     |     |     |     |     |     |  x  |     |     |     |     |     |     |
|  **mik**  |     |     |     |     |     |     |     |     |     |  x  |     |     |  x  |
|  **ger**  |     |     |     |     |     |     |     |     |     |     |     |  x  |     |


## Detailed contributions

| Project task                      | Bot | Developers |
|-----------------------------------|-----|------------|
|Bots integration                   | all |vic         |
|Container bot                      | cnt |vic         |
|2D Delaunay                        | geo |vic         |
|2D Alpha-complex                   | geo |vic         |
|2D meshes                          | geo |vic         |
|Static elasticity FEM assembler    | pde |vic, mvf    |
|Approximated minimum degree        | grp |vic         |
|Spectral bisection                 | grp |vic         |
|Rasterizer                         | g2d |vic         |
|System macros                      | all |mvf         |
|Par parse solvers and structures   | sol |mvf, vic    |
|Domain segmentation                | grp |mvf         |
|Nested dissection                  | grp |mvf         |
|3D Tetrahedral mesher (octree)     | geo |jol         |

## Contributions from 'Public domain code'
- (img) Image reader and writers (PNG, JPG, BMP).
  Sean T. Barrett
  ([GitHub repo](https://github.com/nothings/stb))
- (geo) Robust predicates.
  Jonathan Richard Shewchuk
  ([Web page](https://www.cs.cmu.edu/~quake/robust.html))
    