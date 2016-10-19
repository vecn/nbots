# NBOTS

## Abstract

NBOTS is a scientific computing library focused on parallel computations, fast-prototyping and research.
It has been used in a wide range of architectures and operative systems, such as Android, Mac OS, Linux, Windows and Linux-based supercomputing clusters. 
NBOTS is the short name for 'Numerical Bots', since the library is formed by several interdependent specialized modules, which are referred as bots.

## Compilation

Run `gradle assemble` or use the shell script `./gradlew assemble` (`gradlew.bat` for Windows) if gradle is not installed in your system.

## Numerical bots

The project nbots is a self-contained ecosystem of numerical bots.
A numerical bot is an interdependent module which perform several related numerical tasks.


### Fast reference

|          Name           |                Caption                 |
|-------------------------|----------------------------------------|
| Memory bot              | Handles memory allocations             |
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
### Container bot
### Math bot
### Graphics bot
### Graph bot
### Solver bot
### Optimization bot*
### Statistics bot
### Geometric bot
### Image bot
### Metaheuristics bot*
### PDE bot
### Topopt bot*

## Environment variables
|   Var name  | Bot |    Description      |
|-------------|-----|---------------------|
|NB_FONTS_DIR | g2d | True-type fonts dir |

## Acknowledges

The Center for Research in Mathematics, CIMAT, in many ways has been the backbone of this project, especially **the people** forming the department of Computational Sciences and its post-graduate programs, a big **thank you** for all of them.

## Interesting docs
- [Frequently Asked Questions (FAQ)](FAQ.md)
- [License](LICENSE.md)
- [README for developers](README_DEVELOPERS.md)
- [About developers and contributions](CONTRIBUTIONS.md)