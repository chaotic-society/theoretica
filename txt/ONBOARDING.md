# On-boarding
Welcome to Theoretica! 👋

This quick guide will get you set up and ready for your first contribution. Theoretica is a C++ math library for scientific computing, designed to be simple, elegant, and powerful. It provides numerical methods for calculus, linear algebra, optimization, statistics, and more, all in a header-only library with no dependencies.

## Quick Setup

### Prerequisites
Before diving in, make sure you have:

- A [GitHub](https://github.com) account
- A C++ compiler supporting C++14 or later (e.g. GCC, Clang, MSVC)
- [Git](https://git-scm.com/downloads) or [GitHub Desktop](https://desktop.github.com/download/)
- A text editor or IDE (VS Code, Emacs, Vim, or any you prefer)

### Get the Code
Use Git to download the library code to your local machine:

```bash
git clone https://github.com/chaotic-society/theoretica.git
cd theoretica
make all
```

You should see compilation output followed by test results. This builds all example programs and runs the test suite. If you see errors, make sure your compiler supports C++14 and is properly installed.

If you prefer using **CMake**, that works too:

```bash
cd build
cmake ..
make all
```

If tests pass, you're ready to go! 🎉

## Project Structure

Looking at the repository, you'll see several directories.  Here's what each contains:

```
theoretica/
├── src/                    # Library headers (organized by module)
│   ├── algebra/            # Vectors, matrices, linear algebra
│   ├── autodiff/           # Automatic differentiation
│   ├── calculus/           # Derivatives, integrals, ODEs
│   ├── core/               # Elementary functions (sqrt, exp, sin...)
│   ├── optimization/       # Root-finding, minimization
│   ├── statistics/         # Statistical functions
│   └── theoretica.h        # Main include file
│
├── test/                   # Tests and benchmarks
│   ├── prec/               # Precision tests
│   └── benchmark/          # Performance benchmarks
│
├── examples/               # Example programs
├── txt/                    # Guides and documentation
└── build/                  # Build configuration (CMake)
```

**Key points:**
- **`src/` is where the library code lives**, organized by functionality
- **`test/` contains tests** that verify the library works correctly
- **`examples/` has programs** demonstrating how to use various features
- **`txt/` holds documentation** like coding standards and guides

## Your First Contribution
The easiest ways to start contributing: 

### 1. Improve Documentation
Browse `src/` and look for functions with minimal comments.  When you find one: 
- Read the code to understand what it does
- Add or improve the Doxygen documentation

See the [Contributing Guide](https://github.com/chaotic-society/theoretica/blob/master/txt/CONTRIBUTING.md#contributing-with-documentation) for documentation format.

### 2. Write Tests
Look through `test/prec/` to see existing tests and `src/` to see library code, then: 
- Find functions lacking test coverage
- Add tests for edge cases
- Test mathematical identities

See the [Testing Guide](https://github.com/chaotic-society/theoretica/blob/master/txt/TESTING.md) for details.

### 3. Pick a "Good First Issue"
Visit the [Issues page](https://github.com/chaotic-society/theoretica/issues) and filter by `good first issue`. These are designed for newcomers and include context to help you get started.

For the complete contribution process, see the [Contributing Guide](https://github.com/chaotic-society/theoretica/blob/master/txt/CONTRIBUTING.md).

## Module Organization
Each subdirectory in `src/` is a module focused on a specific area:

- **algebra**: Vectors, matrices, linear systems, decompositions
- **calculus**: Derivatives, integrals, differential equations
- **optimization**: Finding roots and minima of functions
- **statistics**: Mean, variance, regression, distributions

When you're looking for where to add a feature or fix a bug, think about which module it belongs to.

## Getting Help

- **Ask on issues**: Comment directly on the issue you're working on
- **Read the guides** in `txt/` for detailed information
- **Check existing code**: Look at similar functions for guidance

Don't hesitate to ask questions, we're here to help!

## Next Steps

Ready to dive deeper? 

- **[Contributing Guide](https://github.com/chaotic-society/theoretica/blob/master/txt/CONTRIBUTING.md)**: Complete guide to all contribution types
- **[Testing Guide](https://github.com/chaotic-society/theoretica/blob/master/txt/TESTING. md)**: Learn precision testing with Chebyshev
- **[Coding Standard](https://github.com/chaotic-society/theoretica/blob/master/txt/CODING_STANDARD.md)**: Code style conventions
- **[API Documentation](https://chaotic-society.github.io/theoretica)**: Complete library reference

Welcome aboard! We're excited to see your contributions. 🚀
