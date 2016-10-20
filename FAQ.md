# Frequently Asked Questions

## About the project

### Which platforms are supported?
- Linux
- Mac OS
- Windows
- Android

### How to compile the library?
In the directory of the project run the script

`./gradlew assemble`

(Use `gradlew.bat` if you are on windows).

### Why C is the programming language?

The main reasons are:

- C can interact with any programming language using a binding,
  but not in the other way around.
- Every computer architecture has at least a minimal C compiler.
- Every new technology is C-compatible.
- Standards for High Performance Computing are integrated directly in C,
  such as MPI, CUDA-C, OpenMP and Posix-threads, but also the most used
  scientific-computing frameworks, such as OpenGL and Frama-C.
- Easily coupled to embedded systems.
- As any native language, C allows a fine tuning for memory allocation
  and thread execution.
    - **Why C instead of C++?**
      Mainly because every C program can be compiled for any C++ compiler,
      but not in the other way around.
      About C++, we miss the use of namespaces and the exceptions handling.
    - **Why C instead of assembler?**
      Because we are focused on high level numerical procedures.

> "With great power comes great responsibility", Uncle Ben

## Memory bot

## Container bot

## Math bot

## Graphics bot

## Graph bot

## Solver bot

## Optimization bot

## Statistics bot

## Geometric Bot

### Why does the model can not be meshed?
 Use the function `nb_model_verify_consistence()` to check the
 model.

## Image bot

## Metaheuristics bot

## PDE Bot

### What numerical method are used to solve PDE?

- Finite Element Method.
- Control Volume Function Approximation (Finite Volume).
- Discontinuous Galerkin FEM.

### What kind of differential equations can be solved?
Current implementations:

- Stress analysis (solid Mechanics, elasticity).

Coming soon:

- Laplace/Poisson equation (heat analysis, diffusive models).
- Maxwell equations (electromagnetic analysis).

### How to define and run customized FEM simulations?
Please choose a tutorial:

- [Stress analysis](tutorials/stress_analysis.md)
- Heat analysis
- Electromagnetic analysis

## Topopt bot