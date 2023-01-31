# Coding standard
If you contribute to the project, please follow this coding standard to help improve code quality, readability and uniformity.

## Scopes
- Declare all functions and classes inside the `theoretica` namespace and eventually other sub-namespaces (only if needed)
- Do not declare global variable unless truly necessary. Only global constants are accepted and should be declared with capital letters only.

## Naming conventions
- Functions and classes should be named using **snake case** (e.g. `function_name`)
- The name of a function should be specific enough to explain what it does but not too long. Simple shorthands may be used to shorten names if their meaning is evident (e.g. "gen" for "generate").
- The name of a class should refer to the mathematical concept it encapsulates
- Constants should be named with all capital letters and underscores for spaces (e.g. `CONSTANT_NAME`)
- Templated types should be named using **camel case** (e.g. `template<typename ThisIsAType>`)

## Indentation
- Use tabs instead of spaces
- Add a space after every comma (e.g. ", ")

## Documentation
Functions and classes should be documented using Doxygen documentation format, using @param and @return to explain what the arguments are and what the function returns.
For example:
```cpp
/// Do absolutely nothing
/// @param x This is useless
/// @return Literally zero
inline int nothing(int x) {
  return 0;
}
```
Descriptions should explain how the function behaves, the underlying algorithm and error states or exceptions.

## Test cases
- Write test cases for the code you write or at least leave information on how it could be tested

## Other considerations
- All functions except class constructors should be declared `inline` if implemented inside headers
- Do not use the GOTO statement
- Procedures (deterministic functions without side effects, like in functional languages) should be preferred over functions whenever possible for mathematical algorithms and methods
