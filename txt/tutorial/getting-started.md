@page tut_getting_started Getting Started

This tutorial shows a minimal setup and first program using Theoretica.

## Prerequisites
To use Theoretica, you're going to need the same tools that you would use to program in C++:

- A C++ compiler, supporting at least C++14
- Access to the Theoretica headers in your include path

### Getting the library
The easiest way to get Theoretica's source files is to download a ZIP archive of the code [from the repo](https://github.com/chaotic-society/theoretica/archive/refs/heads/master.zip) and extract it in a convenient place, which can also be your project's folder. If you have a working Git installation, you can also clone the repo instead of downloading an archive:

```bash
git clone https://github.com/chaotic-society/theoretica.git
```

### Including the library
The next step is to include the library's header files inside your project. To make sure your compiler finds the header files, you generally have to add an argument to your compiler or configure your IDE. On GCC, you can use the following argument:

```bash
-I/path/to/theoretica/src
```

while on Visual C++ you can use:

```bash
/I"C:\path\to\theoretica\src"
```

If you're using your IDE to compile your programs, you should look into the **build configuration** of your project. Note that if you placed Theoretica's source code inside your project's folder, you should be able to already include it with a relative path (e.g. `#include "../theoretica/src/theoretica.h"`), although an installation in another folder may be cleaner.

### Check your installation
Once you can include the header files in your project, try to run the following example code to see if it compiles without errors:

```c++
#include "theoretica.h"
using namespace th;

int main() {
    
    vec3 v = {1, 2, 3};

    io::println(v);
}
```

## Using specific modules
While the generic header `theoretica.h` can be useful, it can also slow down compilation time, especially if you only use a limited set of functionalities. If this is the case, you can simply include the specific modules that you need, for example `#include "calculus/ode.h"` if you need to solve ODEs.

## OpenMP support
Theoretica automatically detects OpenMP when it is enabled, but you may want to enable it manually in your build configuration if it's not already on. OpenMP speeds up computation by parallelizing code, using automatic dispatch. If OpenMP interferes with your program, you can define `DISABLE_OPENMP` to disable it in Theoretica.

## HDF5 support.
If you need to work on HDF5 files using Theoretica (with `io/hdf5.h`), you need to enable it in your build system or configuration. This is done by including and linking the HDF5 C library, which can be installed for example through `vcpkg`. HDF5 can be a heavy dependency and is recommended only if you plan on using this file format, but is absolutely not needed to use the rest of the library.
