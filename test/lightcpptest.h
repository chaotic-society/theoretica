/*
MIT License

Copyright (c) 2017 Mattia Isgrò

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef LIGHTCPPTEST_H
#define LIGHTCPPTEST_H
#include <iostream>


#define TEST_STARTUP()	unsigned int final_result = 0;	\
						std::cout << "Starting testing with LightCppTest..." << std::endl


#define TEST_BEGIN_MODULE(_module_name)	{										\
											auto module_name = #_module_name;	\
											std::cout << "Testing module " <<	\
											module_name << "...\n" << std::endl

#define TEST_END_MODULE()	std::cout << "Ending test on " << module_name << " module...\n" << std::endl;	\
							}


#define TEST_EXIT()	std::cout << "Finished testing: " << final_result << " tests failed" << std::endl;	\
					return final_result


template<typename ReturnType>
struct Function {
	typedef ReturnType Type;
};

template<typename ReturnType, typename ...Parameters>
inline Function<ReturnType> _extract_retype(ReturnType(*func_ptr)(Parameters...)) {
	return Function<ReturnType>();
}
#define FUNCTION_RETURN_TYPE(function_name) decltype(_extract_retype(&function_name))::Type


#define TEST_BEGIN(_function_name)	{																		\
										auto function_name = #_function_name;								\
										auto function_pointer = &_function_name;							\
										FUNCTION_RETURN_TYPE(_function_name) result;							\
										std::cout << "\tTesting " << function_name << "..." << std::endl;	\
										unsigned int function_final_result = 0


#define TEST_BEGIN_VOID(_function_name)	{																	\
										auto function_name = #_function_name;								\
										auto function_pointer = &_function_name;							\
										std::cout << "\tTesting " << function_name << "..." << std::endl;	\
										unsigned int function_final_result = 0


#define TEST_BEGIN_OVERLOAD(_function_name, return_type, ...)	{											\
										auto function_name = #_function_name;								\
										return_type(*function_pointer)(__VA_ARGS__) = &_function_name;		\
										return_type result;													\
										std::cout << "\tTesting " << function_name <<						\
										"(" << #__VA_ARGS__ << ")" <<										\
										 "..." << std::endl;												\
										unsigned int function_final_result = 0


#define TEST_BEGIN_OVERLOAD_VOID(_function_name, return_type, ...)	{										\
										auto function_name = #_function_name;								\
										return_type(*function_pointer)(__VA_ARGS__) = &_function_name;		\
										std::cout << "\tTesting " << function_name <<						\
										"(" << #__VA_ARGS__ << ")" <<										\
										 "..." << std::endl;												\
										unsigned int function_final_result = 0


#define TEST_END()		if(!function_final_result) std::cout << "\tTests on " << function_name << " succeded\n";	\
						else std::cout << "\t" << function_final_result <<											\
							" tests on " << function_name << " failed\n";											\
						std::cout << "\tEnding test on " << function_name << "...\n" << std::endl;					\
					}


#define TEST_EXEC(...)	result = function_pointer(__VA_ARGS__)


#define TEST_EXEC_VOID(...)	function_pointer(__VA_ARGS__)


#define TEST_EQUALS(expected)	if(result != expected)	{				\
									final_result++;						\
									function_final_result++;			\
									std::cout << "\tTest failed:\n" <<	\
									"\t\tExpected value: " << expected << "\n\t\tFunction returned: " <<	\
									result << "\n" << std::endl;	\
								}


#define TEST_DISEQUALS(unexpected)	if(result == unexpected)	{			\
										final_result++;						\
										function_final_result++;			\
										std::cout << "\tTest failed:\n" <<	\
										"\t\tFunction returned an unexpected value: " << result << "\n" << std::endl;	\
									}


#define TEST_MANUAL_EQUALS(value, expected)	if(value != expected)	{				\
												final_result++;						\
												function_final_result++;			\
												std::cout << "\tTest failed:\n" <<	\
												"\t\tExpected value: " << expected <<	\
												"\n\t\tFunction returned: "			\
												<< value << "\n" << std::endl;		\
											}


#define TEST_MANUAL_DISEQUALS(value, unexpected)	if(value == unexpected)	{			\
													final_result++;						\
													function_final_result++;			\
													std::cout << "\tTest failed:\n" <<	\
													"\t\tFunction returned an unexpected value: "	\
													<< value << "\n" << std::endl;					\
												}


#define TEST_ERROR(description) { std::cout << "\tTest failed:\n\t\t" << std::string(description) << "\n" << std::endl;	\
								final_result++;				\
								function_final_result++; }


#endif
