# Frequently Asked Questions

## About the project

### Why C?

The main reasons are:
- C can interact with any programming language using a binding,
  but not in the other way around.
- Every computer architecture has at least a minimal C compiler.
- Every new technology is C-compatible.
- Standards for High Performance Computing are integrated directly in C,
  such as MPI, CUDA-C, OpenMP and Posix-threads, but also the most used
  scientific-computing frameworks, such as OpenGL and Frama-C.
- As any native language, C allows a fine tuning for memory allocation
  and thread execution.
  ("With great power comes great responsibility", Uncle Ben)

But the truth is that there is not a correct answer for this question because
every programmer has his own favorite language.

** Why C instead of C++?**
Mainly because every C program can be compiled for any C++ compiler,
but not in the other way around.
About C++, we miss the use of namespaces and the exceptions handling.

** Why C instead of assembler?**
Because we are not insane programmers.

## Geometric Bot

### Why does the model can not be meshed?
 Use the function `nb_model_verify_consistence()` to check the
 model.
