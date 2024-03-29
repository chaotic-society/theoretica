# Contributing Guide
Contributions are welcome and appreciated. You don't necessarily have to contribute with code, as new ideas, algorithms and suggestions can be as useful as new code.

## Code of Conduct
If you participate to this project, you are expected to follow the [Code of Conduct](https://github.com/chaotic-society/theoretica/blob/master/CODE_OF_CONDUCT.md).
Feel free to suggest modifications to it. By participating you also accept the library's [License](https://github.com/chaotic-society/theoretica/blob/master/LICENSE).

## What is to be done?
The first place to start to look for new features to implement or research is the [Issues](https://github.com/chaotic-society/theoretica/issues) page.
There you can find a list of features which need to be implemented or improved, and bugs to be fixed.
The description of an issue gives information about what needs to be done.
Pinned issues are the most important features that need to be implemented.

## Suggesting new features
You can also suggest new features to be added to the library.
You can do so by opening an issue and explaining what is the new features, how it may work and why it should be added.
It might help to add some resources to further read about the subject if it refers to some algorithm, method or mathematical concept.
Please make sure that your new issue is **not a duplicate** before opening it.

## Contributing with research
The library needs active research of algorithms, approximations and new concepts. The [Bibliography](https://github.com/chaotic-society/theoretica/blob/master/docs/txt/BIBLIOGRAPHY.md) contains the literature used in the development of the library. You can contribute greatly by opening a new [issue](https://github.com/chaotic-society/theoretica/issues) about your idea, where you can describe how it might improve
the library. You can also write in the [Discussion](https://github.com/chaotic-society/theoretica/discussions) page if you want to open a discussion about a certain topic and confront your ideas with others.

## Writing new code
If you want to contribute by writing new code, make sure to add it in the right place and to follow the [Coding Standard](https://github.com/chaotic-society/theoretica/blob/master/docs/txt/CODING_STANDARD.md).
Please read the [Documentation](https://chaotic-society.github.io/theoretica) to know about the directory structure of the library.
In the case that the new feature does not have an already existing header file to be added to, you can create a new header file following this template:
```cpp
#ifndef THEORETICA_FILENAME_H
#define THEORETICA_FILENAME_H

namespace theoretica {
  
  // Your code here...
  
}

#endif
```

## Writing test cases
Test cases are just as important as new code, because they ensure that the library works correctly! You can help writing test cases by adding new cases to existing test programs, in the folder `test`, or by writing entirely new test programs for different parts of the library. The standard way is to create a different program for each module of the library (a "module" being symbolized by a different folder in the `src` directory, e.g. `src/complex` is the `complex` module). You can follow this template to start writing test cases:

```cpp
#include "theoretica.h"
#include "chebyshev/prec.h"

using namespace chebyshev;
using namespace theoretica;


prec::estimate_result test_custom(interval k, Real tol, unsigned int n) {
	// Custom test
}


int main(int argc, char const *argv[]) {

	prec::state.outputFolder = "test/";
	
	prec::setup("module name");

		prec::estimate(...);

		prec::equals(...);

	prec::terminate();
}
```

### Guidelines for code
- Follow the [Coding Standard](https://github.com/chaotic-society/theoretica/blob/master/docs/txt/CODING_STANDARD.md).
- Make sure to run `make all` before making a pull request, to ensure that your commit does not break other code.
- Write test cases for your code or at least leave information on how it could be tested and its expected behaviour.
- Comment your code to make clear what it does, eventually adding links to resources.

## Pull request
When you are sure about your contribution, make a pull request preferably following the [Pull Request Template](https://github.com/chaotic-society/theoretica/blob/master/.github/PULL_REQUEST_TEMPLATE.md) so that it may be reviewed.
**Thank you for your contribution!**
