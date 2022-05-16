# Contributing Guide
Contributions are welcome and appreciated. You don't necessarily have to contribute with code, as new ideas, algorithms and suggestions can be as useful as new code.

## Code of Conduct
If you participate to this project, you are expected to follow the [Code of Conduct](https://github.com/chaotic-society/uroboro/blob/master/CODE_OF_CONDUCT.md).
Feel free to suggest modifications to it. By participating you also accept the library's [License](https://github.com/chaotic-society/uroboro/blob/master/LICENSE).

## What is to be done?
The first place to start to look for new features to implement or research is the [Issues](https://github.com/chaotic-society/uroboro/issues) page.
There you can find a list of features which need to be implemented or improved, and bugs to be fixed.
The description of an issue gives information about what needs to be done.
Pinned issues are the most important features that need to be implemented.

## Suggesting new features
You can also suggest new features to be added to the library.
You can do so by opening an issue and explaining what is the new features, how it may work and why it should be added.
It might help to add some resources to further read about the subject if it refers to some algorithm, method or mathematical concept.
Please make sure that your new issue is **not a duplicate** before opening it.

## Contributing with research
The library needs active research of algorithms, approximations and new concepts.
You can contribute greatly by opening a new [issue](https://github.com/chaotic-society/uroboro/issues) about your idea, where you can describe how it might improve
the library. You can also write in the [Discussion](https://github.com/chaotic-society/uroboro/discussions) if you want to open a discussion about a certain topic and confront your ideas with others.

## Writing new code
If you want to contribute by writing new code, make sure to add it in the right place.
Please read the [Documentation](https://chaotic-society.github.io/uroboro) to know about the directory structure of the library.
In the case that the new feature does not have an already existing header file to be added to, you can create a new header file following this template:
```cpp
#ifndef UROBORO_FILENAME_H
#define UROBORO_FILENAME_H

namespace uroboro {
  
  // Your code here...
  
}

#endif
```

### Guidelines for code
- Use **snake case** (e.g. `weighted_mean`, not `weightedMean`)
- Use **tabs** instead of spaces.
- Make sure to run `make all` before making a pull request, to ensure that your commit does not break other code.
- Write test cases for your code or at least leave information on how it could be tested and its expected behaviour.
- Comment your code to make clear what it does, eventually adding links to resources.

## Pull request
When you are sure about your contribution, make a pull request preferably following the [Pull Request Template](https://github.com/chaotic-society/uroboro/blob/master/.github/PULL_REQUEST_TEMPLATE.md) so that it may be reviewed.
**Thank you for your contribution!**
