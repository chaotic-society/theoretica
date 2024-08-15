# Contributing Guide
Contributions are welcome and appreciated. You don't necessarily have to contribute with code, as new ideas and perspectives can be as useful as code contributions. The main ways to contribute are by writing code, tests and documentation or by researching algorithms. The following sections introduce how to do so.

### Code of Conduct
If you participate to this project, you are expected to follow the [Code of Conduct](https://github.com/chaotic-society/theoretica/blob/master/CODE_OF_CONDUCT.md).
Feel free to suggest modifications to it. By participating you also accept the library's [License](https://github.com/chaotic-society/theoretica/blob/master/LICENSE).

## What is to be done?
The best place to look for new features to implement or research is the [Issues](https://github.com/chaotic-society/theoretica/issues) page.
There, you can find a list of features which need to be implemented or improved, and bugs to be fixed. The description of an issue gives information about the feature and how it may be implemented. Pinned issues contain the most important tasks that need addressing and are more urgent. 

### Suggesting new features
You can freely suggest new features, doing so by opening a new issue describing your idea and how it may be implemented, how it will improve the project in your opinion and any useful resources. Please make sure to check that your new issue is not a duplicate before opening a new one.

## Contributing with code

Writing new code or improving the existing one is a great way to push the project forward! To contribute by programming, please follow these steps:
- Decide which issue to implement (or open a new one)
- Write your new implementation or changes to existing code
- Whenever possible, write test cases and documentation for your code
- Make sure your code works by running `make all`
- Open a **pull request** (you can follow the [Pull Request Template](https://github.com/chaotic-society/theoretica/blob/master/.github/PULL_REQUEST_TEMPLATE.md))

Your code will then be reviewed for merging into the codebase. Depending on the task, you may need some technical resources on the feature to implement. If the issue does not include it, feel free to ask! Remember to follow the [Coding Standard](https://github.com/chaotic-society/theoretica/blob/master/CODING_STANDARD.md) while writing your code. If you need to create a new header, you can follow this template:

```cpp
#ifndef THEORETICA_FILENAME_H
#define THEORETICA_FILENAME_H

namespace theoretica {
  
  // Your code here...
  
}

#endif
```
You can also consult the [Documentation](https://chaotic-society.github.io/theoretica) to know more about the code structure and specific functionalities.

## Contributing with research

The library needs continuous research of algorithms, concepts and theories. You may contribute greatly by researching topics listed in the [Issues](https://github.com/chaotic-society/theoretica/issues) page. This generally includes understanding and sharing how to implement a certain algorithm or concept, how it would benefit the project and how it works under the hood. Your expertise may prove fundamental in implementing a certain method. Note that contributing with research requires prior knowledge in the fields of numerical analysis, scientific computing, algorithms and data structures. The [Bibliography](https://github.com/chaotic-society/theoretica/blob/master/BIBLIOGRAPHY.md) page is a good starting point to look for research materials. You can also write in the [Discussion](https://github.com/chaotic-society/theoretica/discussions) page if you want to open a discussion about a certain topic and confront your ideas with others. Also remember to add consulted literature to the bibliography.

## Contributing with tests

Tests are just as important as new code, because they ensure that the library works correctly and highlight bugs right when they pop up. Tests are written using [Chebyshev](https://github.com/chaotic-society/chebyshev), a unit testing framework designed for scientific applications, and are stored in the `test` folder. Tests are organized in different implementation files for each module of the library. For example, the `core` module is tested in `test/test_core.cpp`. When creating a new test program, you can follow this simple template:

```cpp
#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {
	
	prec::setup("example");

		output::state.outputFiles = { "test/prec_example.csv" };

		// prec::estimate(...);

		// prec::equals(...);

	prec::terminate();
}
```

Chebyshev has three different modules `prec`, `benchmark` and `err`, but you will generally only need the first one, to estimate the error on numerical approximations. To know more about writing tests, please have a look at [TESTING.md](https://github.com/chaotic-society/theoretica/blob/master/TESTING.md).

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

Classes may be likewise documented using

```cpp
/// @class MyClass
/// Class description
class MyClass { ... }
```

and namespaces with:

```cpp
/// @namespace GlobalNamespace::Namespace Describe this namespace
namespace GlobalNamespace {
	namespace Namespace {
    	/// ...
    }
}
```

To know more about writing Doxygen documentation, please have a look at this [Doxygen guide](https://www.doxygen.nl/manual/docblocks.html).

### Thank you for your contribution!
