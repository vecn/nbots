# NBOTS

## Abstract

NBOTS is a scientific computing library focused on parallel computations,
fast-prototyping and research.
It has been used in a wide range of architectures and operative systems, such
as Android, Mac OS, Linux, Windows and Linux-based supercomputing clusters. 
NBOTS is the short name for 'Numerical Bots', since the library is formed by
several interdependent specialized modules, which are referred as bots.

## Gallery

<img alt="Fast mesh generation using density maps" src="https://user-images.githubusercontent.com/16961118/27419760-2b62f220-56e8-11e7-8ed8-fc582058b45b.gif" >
<img alt="Artistic mesh for demonstration" width="350" src="https://user-images.githubusercontent.com/16961118/27419764-2b6d998c-56e8-11e7-84e4-6cd47b8c0e3e.png">
<img alt="Artistic mesh created for a FEM dedicated blog" width="350" src="https://user-images.githubusercontent.com/16961118/27419762-2b6b5fdc-56e8-11e7-924b-b3362b27de98.png" >
<img alt="Partition of a mesh into 16 subdomains" width="250" src="https://user-images.githubusercontent.com/16961118/27419761-2b66364c-56e8-11e7-9023-785de3f2ea4b.png" >
<img alt="Brake system analysed using FEM" width="550" src="https://user-images.githubusercontent.com/16961118/27419754-2b537264-56e8-11e7-925b-cc8622f2ff64.png" >
<img alt="Shape functions used in Topology optim" width="200" src="https://user-images.githubusercontent.com/16961118/27419763-2b6b8e6c-56e8-11e7-91ba-b80b4a3671e8.gif" >
<img alt="Derivative of a square FE" width="200" src="https://user-images.githubusercontent.com/16961118/27419765-2b6f8cba-56e8-11e7-8c71-b458442ee36d.gif" >
<img alt="Topology optimization of a cantilever" src="https://user-images.githubusercontent.com/16961118/27419757-2b5702a8-56e8-11e7-9310-9c8d365e5789.gif" >
<img alt="Ultra-brittle fracture simulated with DEM" src="https://user-images.githubusercontent.com/16961118/27419755-2b55e724-56e8-11e7-8ed2-d530ba2b03c5.gif" >
<img alt="Detonation using a small explosive load" src="https://user-images.githubusercontent.com/16961118/27419756-2b55ee4a-56e8-11e7-9cf3-2e61e356e2f4.gif" >
<img alt="Detonation using a medium explosive load" src="https://user-images.githubusercontent.com/16961118/27419758-2b574966-56e8-11e7-9012-da169c716238.gif" >
<img alt="Detonation using a huge explosive load" src="https://user-images.githubusercontent.com/16961118/27419759-2b57e1a0-56e8-11e7-8287-bc3f600207e7.gif" >

## Compilation

Create a **build** directory and inside the directory run

```
build$> cmake [nbots dir]
build$> make
```

where **nbots dir** is the route to the project.

## Numerical bots

The project nbots is a self-contained ecosystem of numerical bots.
A numerical bot is an interdependent module which perform several related
numerical tasks.


### Fast reference

|          Name           |                Caption                 |
|-------------------------|----------------------------------------|
| Memory bot              | Handles memory allocations             |
| IO bot                  | IO structures and functions            |
| Container bot           | Process dynamic collections            |
| Math bot                | Standard numerical procedures          |
| Graphics bot            | 2D-vector drawing functions            |
| Graph bot               | Graph-specific routines                |
| Solver bot              | Linear algebra operations              |
| Optimization bot        | Constrained and unconstrained optim    |
| Statistics bot          | Pseudo-random number generators        |
| Geometric bot           | Implements geometric algorithms        |
| Image bot               | Image processing routines and filters  |
| Metaheuristics bot      | Black-box and combinatorial optim      |
| PDE bot                 | Solve PDE using FEM & CVFA             |
| Topopt bot              | Topology optimization problems         |

### Memory bot

This module is provided for memory management.
All memory allocations within the project must be done using this bot.

### Container bot

The Container bot is a collection of data structures for handling dynamic
collections.
The computational complexity must be considered when designing an algorithm.
The following tables summarizes the memory required per item inserted and
the worst case execution bound per operation (big O notation).

|Container | Implementation | Bytes x Item | Pointers x item | Sorted? | Pointer repetition? |
|-------|----------------------|:---:|:---:|:---:|:---:|
|Queue	| Circular Linked-list |  0  |  2  | no  | yes |
|Stack	| Circular Linked-list |  0  |  2  | no  | yes |
|Sorted	| AVL                  |  4  |  3  | yes | no  |
|Hash	| Array of queues      |  0  |  2  | no  | yes |
|Heap	| Half sorted BT       |  1  |  3  | yes | yes |


|Container | insert() | delete() | exist() | merge() | get_fisrt() | delete_first() | iteration |
|-------|:----:|:----:|:----:|:------:|:----:|:----:|:----:|
|Queue	|   1  |   N  |   N  |    1   |   1  |   1  |   1  |
|Stack	|   1  |   N  |   N  |    1   |   1  |   1  |   1  |
|Sorted	|log(N)|log(N)|log(N)|N log(N)|log(N)|log(N)|log(N)|
|Hash	|   1  |   1  |   1  |    N   |   N  |   N  |   1  |
|Heap	|   1  |   N  |   1  |    1   |   1  |log(N)|log(N)|

### IO bot

Input Output functions.

### Math bot

Module provided for standard numerical procedures, such as numerical quaratures
for polynomial integration.

### Graphics bot

The Graphics bot is a 2D vector graphics library which can be used to generate
drawings and useful visualizations.

### Graph bot

The graph bot is provided for graph processing, such as
- Labeling
- Spectral clusterization
- Perfect matching
- Coloring

### Solver bot

Linear algebra module focused on parallel computations and sparse matrices.

### Optimization bot

Module provided for classical optimization, constrained and unconstrained.

### Statistics bot

Bot provided for random number generators.

### Geometric bot

The geometrical bot is provided for geometric operations, such as
- Meshing.
- Boolean shape operations.
- Delaunay triangulations.
- Voronoi diagrams.
- Point checks.
- Intersection checks.
- Fast geometrical search.

### Image bot

The image bot is intended for image processing.

### Metaheuristics bot

This module is focused in black box and combinatorial optimization.

### PDE bot

Module specialized in solving partial differential equations (PDE) using the
most accurate numerical methods:

- Finite Element Method (FEM)
- Control Volume Function Approximation (CVFA)

Customized schemes for engineering problems:
- Solid Mechanics.
- Electromagnetic simulations.
- Heat transfer.

### Topopt bot

This bot is provided for "Topology optimizations" problems.

## Environment variables
|   Var name  | Bot |    Description      |
|-------------|-----|---------------------|
|NB_FONTS_DIR | g2d | True-type fonts dir |

## Acknowledges

The Center for Research in Mathematics, CIMAT, in many ways has been the
backbone of this project, especially **the people** forming the department of
Computational Sciences and its post-graduate programs, a big **thank you** for
all of them.

## Interesting docs
- [Frequently Asked Questions (FAQ)](FAQ.md)
- [License](LICENSE.md)
- [README for developers](README_DEVELOPERS.md)
- [About developers and contributions](CONTRIBUTIONS.md)