# Dear developers

I wrote this 'Quick-start guide' as a friendly introduction to development conventions, coding style and other caveats of this library. I hope you find it useful.

## Why C?

The main reasons are:
- C can interact with other programming languages using a binding, but not in the other way around.
- Every computer architecture has at least a minimal C compiler.
- Every new technology is C-compatible.
- Standards for High Performance Computing are integrated directly in C, such as MPI, CUDA-C, OpenMP and Posix-threads, but also the most used scientific-computing frameworks, such as OpenGL.
- C can be compiled for any C++ compiler, but not in the other way around.
- C allows a fine tuning for memory allocation and thread execution.
  ("With great power comes great responsibility", Uncle Ben)

## Coding style

Every one has his own programming style and that is fine, but when working in a sustainable project, it is desirable to follow some coding style and conventions, in order to make the project readable and maintainable for everyone. We follow most of the rules used by the linux kernel developers:

https://www.kernel.org/doc/Documentation/CodingStyle

## Fast checklist

> Produce files with no more than 500 rows (aims a modular code).
> Declare static those functions used just in a single file.
> Add the prefix 'nb_' to those functions visible to every file.
> Write your code using 80 columns at most (increase readability).
> Write functions of 40 rows at most (to compile with the mind).
> 3 identation levels are enough, 4 levels are some times required
  (more than 4 provokes brain damage).
> Use 3 control statements at most per function (for, if, while).
> Make just one memory allocation per function when required.
> Never use malloc(), calloc() and realloc() directly.
  Use those declared on nb/memory_bot/allocate_mem.h
> If the allocated memory will be freed on the same function,
  use the "soft_allocate" version.
> Use 'nb_membank_t' if it is possible.
> Deliver at least one unit test for every function declared in the header 
  (Using CUNIT).

## Library standard
The library uses objects (structs) to handle the data, every object must have at least the following functions:

> nb_[object_name]_get_memsize();
> nb_[object_name]_init(void *ptr_to_object);
> nb_[object_name]_finish(void *ptr_to_object);
> nb_[object_name]_copy(void *ptr_to_object, const void *src);
> nb_[object_name]_clear(void *ptr_to_object);

The create(), destroy() and clone() functions can be built upon these standard functions, but are not included because they are memory exhaustive and decrease the performance on mobile devices.

## Building the library

To build the project we use 'Gradle', one of the most used building tools for big projects.

To build the library run **gradle assemble**.
To build the library and execute the unit tests run **gradle build**.
If you do not have **gradle** installed in your system, use the wrapper, **gradlew.sh** on unix/like systems or **gradlew.bat** on windows.

## Final comments
- This library has been tested on Windows, Max OS, Linux and Android.
- Please read the book 'Clean code' written by Robert Cecil Martin.

Best regards.

Victor Eduardo Cardoso Nungaray
Guanajuato, Gto. Mexico, October 2016.