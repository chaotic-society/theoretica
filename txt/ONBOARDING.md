# On-boarding
This tutorial will guide you through your first steps to contribute to Theoretica,
be it with code, documentation or tests. In any case, your contributions will help shape and improve the project.
First of all, you will need the following things:
- A [Github](https://github.com) account
- A C++ compiler (not needed for documentation writing)
- A [Git](https://git-scm.com/downloads) installation (or [Github Desktop](https://desktop.github.com/download/))

If you haven't, a good idea is to read the [README](https://github.com/chaotic-society/theoretica/blob/master/README.md)
of the library to get a first idea of what's going on.

## First steps
Once you have a working setup, you need to _clone_ the repository of the project to your computer.
Open a terminal window in the folder where you want to copy the code and launch the following command:

```
git clone https://github.com/chaotic-society/theoretica.git
```

This will copy the code to a folder called `theoretica` and you will be able to modify the code and contribute.
Now, to check that everything is working correctly, you can use `make` (or if you prefer, `CMake`) to build the
example programs and tests. Run the following commands:

```
cd theoretica
make all
```

You should be greeted by a long wall of text containing all sorts of function names, those are the test programs of the library!
If you get any errors and are unable to compile the project, copy the output and contact a maintainer to understand what's not working.

## Where should I look?
Looking into the folder, you can see all sorts of things. It contains not only code but also documentation and configuration files.
Let's have a closer look at the project structure:

- The `build` folder contains configuration files and you won't generally need to modify it.

- The `examples` folder contains example programs using the library. These can be helpful to learn how to use it.

- The `src` folder is where the magic is at. It contains all header files of the library and this is where you can modify or extend the codebase.

- The `test` folder contains precision testing programs and benchmarks. This is where you will be able to write tests, in particular the `prec` subdirectory.

- The `txt` folder contains many useful files (just like this one) which will tell you more about different parts of the library.

The documentation is stored directly alongside code, that's why there is no separate folder for writing documentation (the documentation is then automatically generated and uploaded to the `gh-pages` branch).

## How to Write Documentation
> Writing and maintaining documentation is one of the easiest and yet most important things to do in any open source project.
This is a good task to get going with a project as it makes you read code and try to understand it.

To contribute with documentation writing, look through the library code in `src` and find functions which are lacking documentation or have little to no explanation.
Read and try to understand what the function does and how the parameters work (don't worry if you struggle, you can always ask a maintainer),
and then fill in the missing documentation with your new information, like in this example:

```cpp
/// Compute the square root.
real sqrt(real x) {
  ...
}
```

This function has a Doxygen comment block marked by `///` (3 slashes instead of 2), but has only a short description of what it does.
We can add important information using the `@param` (for function arguments) and `@return` (for return values) Doxygen directives:

```cpp
/// Compute the square root of a non-negative real number,
/// using Newton's method or x86 operations if available.
///
/// @param x A non-negative real number
/// @return The square root of x
real sqrt(real x) {
  ...
}
```

In addition, you can use `@tparam` to document template parameters, with a syntax analogous to `@param`.
You can also document classes with `@class` and files with `@file`. Every source file and class should be
documented, explaining their role and functioning.
