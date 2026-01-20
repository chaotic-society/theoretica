# Contributing Guide 🏗️
Welcome to Theoretica! We're happy that you're interested in contributing. Whether you're here to fix a bug, add a feature, improve documentation, or simply learn, you're in the right place.

Contributions come in many forms, and you don't necessarily need to write code to make a meaningful impact. New ideas, perspectives, bug reports, and documentation improvements are just as valuable as code contributions.

## Ways to Contribute
There are several ways you can help improve Theoretica: 

| Contribution Type | Description | Good For |
|-------------------|-------------|----------|
| 🐛 **Bug Reports** | Report issues you've found | Everyone |
| 💡 **Feature Requests** | Suggest new functionality | Everyone |
| 💻 **Code** | Implement features or fix bugs | Developers |
| 🧪 **Tests** | Write tests to verify correctness | Developers & Beginners |
| 📚 **Documentation** | Improve docs and examples | Developers & Beginners |
| 🔬 **Research** | Research algorithms and theory | Mathematicians & Scientists |

### What is to be done?
The best place to find tasks is the [Issues](https://github.com/chaotic-society/theoretica/issues) page. Issues are labeled to help you find what you're looking for: 

- **`good first issue`** — Great for newcomers to the project
- **`help wanted`** — We'd especially appreciate help here
- **`bug`** — Something isn't working correctly
- **`enhancement`** — New features or improvements
- **`documentation`** — Improvements to docs

Each issue description provides context about what needs to be done and may include hints on implementation. Feel free to ask questions in the issue comments if anything is unclear!

### Reporting Bugs
If you have found unexpected behaviors, errors or other inconsistencies while using Theoretica, please report it by creating a new [issue](https://github.com/chaotic-society/theoretica/issues), writing your findings and including code and information about your environment. By doing so, you're helping keep Theoretica bug-free.

### Suggesting New Features
Have an idea for something new? We'd love to hear it! You can [suggest new features](https://github.com/chaotic-society/theoretica/issues/new?assignees=&labels=&projects=&template=feature_request.md&title=) by opening a new issue. Please describe: 

- **What** you'd like to see implemented
- **Why** it would be useful (use cases, motivation)
- **How** it might work (if you have ideas)
- **References** to relevant papers, algorithms, or existing implementations

Don't worry if you're not sure about all the details, we can discuss and refine the idea together!

## Contributing with Code
Writing new code or improving existing implementations is a wonderful way to push the project forward. Here's how to get started: 

### Step-by-Step Process

- [ ] **Find or create an issue**
   - Pick an existing [issue](https://github.com/chaotic-society/theoretica/issues) that interests you, or open a new one
   - Comment on the issue to let others know you're working on it

- [ ] **Fork the repository**
   - [Fork](https://github.com/chaotic-society/theoretica/fork) the repo on GitHub to create your own copy

- [ ] **Clone your fork locally**
   ```bash
   git clone https://github.com/YOUR_USERNAME/theoretica.git
   cd theoretica
   ```
   Alternatively, you can use a Git client such as [GitHub Desktop](https://github.com/apps/desktop) to avoid using terminal commands. Just follow the same steps, but inside the graphical interface.

- [ ] **Create a feature branch (Optional)**
   ```bash
   git checkout -b feature/your-feature-name
   ```

- [ ] **Make your changes**
   - Write your implementation
   - Follow the [Coding Standard](https://github.com/chaotic-society/theoretica/blob/master/txt/CODING_STANDARD.md)
   - Add documentation (see [Doxygen format](#documenting-your-code))
   - Write tests when possible (see [Testing Guide](https://github.com/chaotic-society/theoretica/blob/master/txt/TESTING.md))

- [ ] **Test your changes**
   ```bash
   make all		# Build everything and run tests
   ```

- [ ] **Commit your changes**
   ```bash
   git add .
   git commit -m "Add brief description of your changes"
   ```

- [ ] **Push to your fork**
   ```bash
   git push origin feature/your-feature-name
   ```

- [ ] **Open a Pull Request**
   - Open a Pull Request on the original [repo](https://github.com/chaotic-society/theoretica/compare)
   - Follow the [Pull Request Template](https://github.com/chaotic-society/theoretica/blob/master/.github/PULL_REQUEST_TEMPLATE.md)

### Creating a New Header File
If you need to create a new header file in `src`, follow this template:

```cpp
///
/// @file filename.h Brief description of this file
///

#ifndef THEORETICA_FILENAME_H
#define THEORETICA_FILENAME_H

namespace theoretica {
  
    // Your code here... 
  
}

#endif
```

**Important tips:**
- Use include guards with the `THEORETICA_` prefix
- Document the file purpose with a `@file` directive
- Keep everything inside the `theoretica` namespace

### Documenting Your Code
All public functions and classes should be documented using Doxygen format: 

```cpp
/// A short description of the function.
///
/// @param x Describe what this function argument is
/// @return Describe what the function returns
///
/// A longer description of the function. For example,
/// which exceptions it may throw and under which conditions.
inline real f(real x) {
    // ...
}
```

You can read more about Doxygen documentation [here](https://www.doxygen.nl/manual/docblocks.html).

## Contributing with Testing
Tests are essential for ensuring the library works correctly and for catching bugs early. Even if you don't write new features, contributing tests is valuable!

### Getting Started with Tests
Tests are written using the [Chebyshev](https://github.com/chaotic-society/chebyshev) framework and live in the `test/prec/` directory. Here's a minimal test file:

```cpp
#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;

int main(int argc, char const *argv[]) {
    
    // Create a test context
    auto ctx = prec::make_context("moduleName", argc, argv);
    ctx.settings.outputFiles = { "test/prec/prec_module.csv" };

    // Test exact values
    ctx.equals("funcName", x_approx, x_exact, options);

    // Test over a range
    ctx.estimate("funcName", f_approx, f_exact, options);
    
    return 0;
}
```

For comprehensive testing guidelines, please see the [Testing Guide](https://github.com/chaotic-society/theoretica/blob/master/txt/TESTING.md).

### What to Test

- **Numerical accuracy**: Does the function produce correct results?
- **Edge cases**: How does it handle zeros, negatives, infinities?
- **Boundary conditions**: What happens at the limits of the input domain?
- **Stability**: Does it maintain precision across many operations?

## Contributing with Research
The library relies on continuous research of algorithms, mathematical concepts, and numerical methods. If you have expertise in mathematics, physics, or scientific computing, your research contributions are invaluable!

You may contribute greatly by researching topics listed in the [Issues](https://github.com/chaotic-society/theoretica/issues) page. This generally includes understanding and sharing how to implement a certain algorithm or concept, how it would benefit the project and how it works under the hood. Your expertise may prove fundamental in implementing a certain method. Note that contributing with research requires prior knowledge in the fields of numerical analysis, scientific computing, algorithms and data structures.

- **Find better algorithms**: Research faster or more accurate approaches
- **Provide references**: Share papers, books, or resources relevant to issues
- **Analyze edge cases**: Help identify numerical pitfalls and stability concerns
- **Verify implementations**: Check that our methods match theoretical descriptions

The [Bibliography](https://github.com/chaotic-society/theoretica/blob/master/txt/BIBLIOGRAPHY.md) page is a good starting point to look for research material. You can also write in the [Discussion](https://github.com/chaotic-society/theoretica/discussions) and [Issues](https://github.com/chaotic-society/theoretica/issues) pages if you want to open a discussion about a certain topic and confront your ideas with others. Also remember to add consulted literature to the bibliography.

## Contributing with documentation
Another way to contribute to the project is by writing documentation. The library is documented by writing Doxygen documentation right alongside the code. This makes it easy to continuously update it and store it. The generated HTML documentation is available online at [Documentation](https://chaotic-society.github.io/theoretica) or for download in the `gh-pages` branch. There is no need to manually build the documentation as on each commit it is automatically generated and uploaded by a Github Actions workflow. You may follow these steps to contribute with documentation:

- Find undocumented/not fully documented code or ask on Github for indications
- Read the code thoroughly and understand what it does and how
- Write or improve the documentation using Doxygen notation
- Make a pull request for your changes

Doxygen is really easy to use and boils down to writing comments like this, right before the construct (in the case of functions):

```cpp
/// Short description
///
/// @param var Description of parameter variable
/// @return Description of return value
///
/// Long description of my function
void myFunction() { ... }
```

Classes may be likewise documented using:

```cpp
/// @class MyClass
/// Class description
class MyClass { ... }
```

and namespaces with:

```cpp
/// @namespace global_namespace::namespace Describe this namespace
namespace global_namespace {
	namespace namespace {
    	/// ...
    }
}
```

To know more about writing Doxygen documentation, please have a look at this [Doxygen guide](https://www.doxygen.nl/manual/docblocks.html).

## Code Review Process

After you submit a Pull Request: 

1. **Automated Checks**: Your code will be built and tested automatically
2. **Review**: A maintainer will review your changes
3. **Feedback**: You may receive suggestions or requests for changes
4. **Approval**: Once everything looks good, your contribution will be accepted!

Don't worry if you receive feedback, it's a normal part of the collaborative process, and we're here to help you succeed.

### Thank you for your contribution!
