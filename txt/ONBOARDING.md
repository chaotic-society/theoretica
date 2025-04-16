# On-boarding
This tutorial will guide you through your first steps to contribute to Theoretica,
be it with code, documentation or tests. First of all, you will need the following things:
- A [Github](https://github.com) account
- A C++ compiler (not needed for documentation writing)
- A [Git](https://git-scm.com/downloads) installation (or [Github Desktop](https://desktop.github.com/download/))

## First steps
Once you have a working setup, you need to _clone_ the repository for the project to your computer.
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
