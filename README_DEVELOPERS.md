# Dear developers

I wrote this 'Quick-start guide' as a friendly introduction to development conventions, coding style and other caveats of this library.
I hope you find it useful.

## Coding style

Every one has his own programming style and that is fine, but when working in a sustainable project, it is desirable to follow some coding style and conventions, in order to make the project readable and maintainable for everyone. We follow most of the rules used by the linux kernel [coding style](https://www.kernel.org/doc/Documentation/CodingStyle)

Since the library implementation begins in 2009, old code does not follow all the conventions and coding style discussed here, but we are normalizing all the library as fast as we can.

## Fast checklist

- About dependencies
    * The code must not depend on any library.
    * The code must be ANSI C (ISO C89). 
        - `for` statement does not allow inline variable declaration.
        - `bool` not supported.
- Coding Style and format
    * Write your code using 80 columns at most (increase readability).
    * Write  functions of 40 rows at most (to compile with the mind).
    * 3 indentation levels are enough, 4 levels are some times required
      (more than 4 provokes brain damage).
    * 8 chars per indentation level
      (improves readability; discussed on the CodingStyle document cited).
    * Use 3 control statements at most per function (for, if, while).
    * Use nouns for classes/structures and verbs for functions.
- About comments
    * All comments must be C-style, that is  `/* comment here */`,
      avoid using C++ single line comments `// like this comment`.
    * The files must not contain any comment at the top.
        + Do not paste the License 
	  (The License is on the **git** repository).
        + Do not sign the files with your name and the date
	  (This info is also in the repository,
	   and the **README** must contain all your
	   contributions, acknowledges and credits).
    * Do not let code-lines commented, delete them if they are not
      required anymore (The repository has the historial of your code,
      you can always go back and recover the deleted lines).
    * Do not make silly comments, e.g. int a = 1; /* Adding one to a */
    * Add a comment only if it is necessary for explaining some weird
      operation or an important note (if you need a lot of comments to
      explain your code, then the code must be rewritten in an elegant
      manner).
- Library structure
    * Produce files with no more than 500 rows (aims a modular code).
    * Declare static those functions used just in a single file.
    * Add the prefix 'nb_' to those functions visible to every file.
    * The directory structure of 'include', 'src' and 'utest' dirs must
      be the same. For each header-file under 'include' dir, it must
      exist at least one file with the same name under 'src' and 'utest'.
- Error handling
    * Return an integer error code  (zero if success).
    * Use `nb_assert()` for check developer exclusive errors.
    * For large verifications using only assert, check `#ifndef  NB_DEBUG`
- Memory management
    * Make just one memory allocation per function when required.
    * Never use `malloc()`, `calloc()` and `realloc()` directly.
      Use those declared on **nb/memory_bot/allocate_mem.h**
    * If the allocated memory will be freed on the same function,
      use the `nb_soft_allocate_mem()` version.
    * Use `nb_membank_t` whether it is possible.
    * Deliver at least one unit test for every function declared in the
      header (Using CUNIT).

## Library standard
The library uses objects (structs) to handle data, every object must have at least the following functions:

- `nb_[object_name]_get_memsize();`
- `nb_[object_name]_init(void *ptr_to_object);`
- `nb_[object_name]_finish(void *ptr_to_object);`
- `nb_[object_name]_copy(void *ptr_to_object, const void *src);`
- `nb_[object_name]_clear(void *ptr_to_object);`

The `create()`, `destroy()` and `clone()` functions can be built upon these standard functions, but are not included because they are memory exhaustive and decrease the performance on mobile devices.

## Building the library

To build the project we use **Gradle**, one of the most used building tools for big projects.

To build the library run `gradle assemble`.
To build the library and execute the unit tests run `gradle build`.
If you do not have **gradle** installed in your system, use the wrapper, `./gradlew.sh` on unix/like systems or `gradlew.bat` on windows.

## Three laws of Test Driven Development (TDD)
1. You are not allowed to write any production code until you have first written a failing unit test.
2. You are not allowed to write more of a unit test than is sufficient to fail (and not compiling is failing).
3.  You are not allowed to write more production code that is sufficient to pass the currently failing unit test.

## Final comments
- This library has been tested on Windows, Max OS, Linux and Android.
- We use CUNIT to perform the unit tests via Gradle.
- We will incorporate Frama-C for code analysis.
- Please read the book 'Clean code' written by Robert Cecil Martin
  (and if possible read also 'Clean coder')

> "Getting software to work and getting software clean are two different activities"
Robert C. Martin

Best regards,

Victor Eduardo Cardoso Nungaray.

Guanajuato, Gto. Mexico, October 18, 2016.