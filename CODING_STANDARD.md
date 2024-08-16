# Coding standard
If you contribute to the project, please follow this coding standard to help improve code quality, readability and uniformity. Thanks for your contribution!

## Scopes
- Declare all functions and classes inside the `theoretica` namespace and eventually other sub-namespaces (only if needed)
- Do not declare global variables unless truly necessary. Only global constants are accepted and should be declared with capital letters only.

## Naming conventions
- Functions and classes should be named using **snake case** (e.g. `function_name`)
- The name of a function should be specific enough to explain what it does but not too long. The naming convention for functions is `action_method_modifier`. For example, `algebra::decompose_lu_inplace` computes the LU matrix decomposition in-place. Abbreviations may be used when the modifier or method is long.
- The name of a class should refer to the mathematical concept, algorithm or idea it encapsulates. For example, the classes `dual` and `complex` implement the corresponding algebra elements.
- Constants should be named with all capital letters and underscores for spaces (e.g. `CONSTANT_NAME`). Configuration constants follow the naming convention `THEORETICA_<CONTEXT_<NAME>`.
- Macros should be used only when necessary, following the naming convention `TH_<NAME>`.
- Templated types should be named using **camel case** (e.g. `template<typename ThisIsAType>`)

## Indentation
- Use tabs instead of spaces
- Add a space after every comma (e.g. ", ")

## Documentation
Functions and classes should be documented using Doxygen documentation format, using `@param` and `@return` to explain what the arguments are and what the function returns.
For example:
```cpp
/// Do absolutely nothing
///
/// @param x This is useless
/// @return Literally zero
///
/// Long description ...
inline int nothing(int x) {
  return 0;
}
```
Descriptions should explain how the function behaves, the underlying algorithm and error states or exceptions. Latex formulas may be written using opening and closing `\f$`, for example `\f$x = \pi\f$`

## Test cases
Write test cases for the code you write or at least leave information on how it could be tested. Tests are written in `test/test_<module>.cpp` for each module, please make sure to include your tests in the right file.

## Other considerations
- All functions except class constructors should be declared `inline` if implemented inside headers
- Do not use the GOTO statement
- Procedures (deterministic functions without side effects, like in functional languages) should be preferred over normal functions whenever possible, for numerical methods
- Do not use the `using namespace` directive inside a header file, as it pollutes the global namespace. It may be used in specific implementation files like test units or benchmarks
