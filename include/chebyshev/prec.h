
///
/// @file prec.h Precision estimation of real functions
///

#pragma once

#include "core/common.h"
#include "core/interval.h"

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <map>
#include <array>
#include <fstream>
#include <string>


#ifndef CHEBYSHEV_INTEGRAL_ITER
#define CHEBYSHEV_INTEGRAL_ITER 10000
#endif


#ifndef CHEBYSHEV_TOLERANCE
#define CHEBYSHEV_TOLERANCE 0.00000001
#endif


namespace chebyshev {
	
	namespace prec {
		
		struct estimate_request;
		struct equation_request;
		struct estimate_result;
		struct equation_result;

		/// Global state of precision testing
		struct prec_state {
			
			/// Name of the module being tested
			std::string moduleName = "unknown";
			
			/// Print to standard output or not
			bool quiet = false;

			/// Write to standard output only failed/not passed estimates?
			bool estimateOnlyFailed = false;

			/// Write to standard output only failed/not passed equations?
			bool equalsOnlyFailed = false;

			/// Output to file?
			bool outputToFile = true;

			/// Output file
			std::ofstream outputFile;

			/// Relative or absolute path to output folder
			std::string outputFolder = "";

			/// Total number of tests run
			uint32_t totalTests = 0;

			/// Number of failed tests
			uint32_t failedTests = 0;

			/// Default number of iterations for integral quadrature
			uint32_t defaultIterations = CHEBYSHEV_INTEGRAL_ITER;

			/// Default tolerance on max absolute error
			Real defaultTolerance = CHEBYSHEV_TOLERANCE;

			/// Recorded estimation requests
			std::vector<estimate_request> estimation_requests;

			/// Recorded equation requests
			std::vector<equation_request> equation_requests;

			/// Results of precision testing
			std::map<std::string, std::vector<estimate_result>> estimation_results;
			
			/// Results of equations
			std::map<std::string, std::vector<equation_result>> equation_results;

		} state;


		/// @class estimate_request A precision estimation request
		struct estimate_request {
			
			/// Uniquely identifying function name
			std::string funcName = "unknown";

			/// The function to estimate
			RealFunction func = nullptr;

			/// A function returning the expected output
			RealFunction funcExpected = nullptr;

			/// Requested estimation intervals
			std::vector<interval> intervals;

			/// Precision testing tolerance on max absolute error
			Real tolerance = state.defaultTolerance;

			/// Number of iterations for integral quadrature
			uint32_t iterations = state.defaultIterations;

			/// Print to standard output or not
			bool quiet = false;
		};


		/// @class equation_request An equation request
		struct equation_request {

			/// Uniquely identifying function name
			std::string funcName = "unknown";

			/// Evaluated value
			Real evaluated;

			/// Expected value
			Real expected;

			/// Tolerance
			Real tolerance = state.defaultTolerance;

			/// Print to standard output or not
			bool quiet = false;
		};


		/// @class estimate_result The result of error estimation
		struct estimate_result {
			
			/// Uniquely identifying name of the function
			std::string funcName = "unknown";

			/// Interval of estimation
			interval k;

			/// Tolerance on the max absolute error
			Real tolerance;

			/// Estimated maximum absolute error on interval
			Real max_err;

			/// Estimated mean error on interval
			Real mean_err;

			/// Estimated RMS error on interval
			Real rms_err;

			/// Estimated relative error on interval
			Real rel_err;

			/// Estimated absolute error on interval
			Real abs_err;

			/// Did the test fail?
			bool failed;

			/// Print to standard output or not
			bool quiet = false;

			/// Total number of iterations for integral quadrature
			uint32_t iterations;
		};


		/// @class equation_result The result of equation checking
		struct equation_result {

			/// Uniquely identifying function name
			std::string funcName = "unknown";

			/// Evaluated value
			Real evaluated;

			/// Expected value
			Real expected;
			
			/// Absolute difference between expected and evaluated values
			Real diff;

			/// Tolerance on the absolute difference
			Real tolerance;

			/// Did the test fail?
			bool failed;

			/// Print to standard output or not
			bool quiet = false;
		};


		/// Register a function for error estimation
		inline void estimate(
			std::string name,
			RealFunction fApprox,
			RealFunction fExp,
			interval k,
			Real tolerance = state.defaultTolerance,
			bool quiet = false,
			unsigned int n = state.defaultIterations) {

			estimate_request r;
			r.funcName = name;
			r.func = fApprox;
			r.funcExpected = fExp;
			r.intervals.push_back(k);
			r.tolerance = tolerance;
			r.quiet = quiet;
			r.iterations = n;

			state.estimation_requests.push_back(r);
		}


		/// Register a function for error estimation on multiple intervals
		inline void estimate(
			std::string name,
			RealFunction fApprox,
			RealFunction fExp,
			std::vector<interval> intervals,
			Real tolerance = state.defaultTolerance,
			bool quiet = false,
			unsigned int n = state.defaultIterations) {

			estimate_request r;
			r.funcName = name;
			r.func = fApprox;
			r.funcExpected = fExp;
			r.intervals = intervals;
			r.tolerance = tolerance;
			r.quiet = quiet;
			r.iterations = n;

			state.estimation_requests.push_back(r);
		}


		/// Estimate the precision of a real function on a single interval
		inline estimate_result compute_estimate(
			std::string name,
			RealFunction fApprox,
			RealFunction fExp,
			interval k,
			Real tolerance = state.defaultTolerance,
			bool quiet = false,
			unsigned int n = state.defaultIterations) {

			estimate_result result;

			Real sum = 0;
			Real sum_sqr = 0;
			Real sum_abs = 0;
			Real max = 0;

			Real measure = k.length();
			Real dx = measure / n;
			Real x, coeff;

			Real diff = std::abs(fApprox(k.a) - fExp(k.a));

			sum += diff;
			sum_sqr += diff * diff;
			sum_abs += std::abs(fExp(k.a));
			max = diff;

			for (unsigned int i = 1; i < n; ++i) {

				x = k.a + i * dx;
				diff = std::abs(fApprox(x) - fExp(x));

				if(diff > max)
					max = diff;

				if(i % 2 == 0)
					coeff = 2;
				else
					coeff = 4;

				sum += coeff * diff;
				sum_sqr += coeff * diff * diff;
				sum_abs += coeff * fExp(x);
			}

			diff = std::abs(fApprox(k.b) - fExp(k.b));

			sum += diff;
			sum_sqr += diff * diff;
			sum_abs += std::abs(fExp(k.b));
			
			if(diff > max)
				max = diff;

			result.funcName = name;
			result.quiet = quiet;
			result.k = k;
			result.iterations = n;
			result.abs_err = sum;
			result.max_err = max;
			result.mean_err = (sum * dx / 3.0) / measure;
			result.rms_err = std::sqrt((sum_sqr * dx / 3.0) / measure);
			result.rel_err = (sum * dx / 3.0) / (sum_abs * dx / 3.0);
			result.tolerance = tolerance;

			if(result.max_err > tolerance) {
				state.failedTests++;
				result.failed = true;
			} else {
				result.failed = false;
			}

			state.estimation_results[result.funcName].push_back(result);
			state.totalTests++;

			return result;
		}


		inline std::vector<estimate_result> compute_estimate(estimate_request r) {

			std::vector<estimate_result> res;

			for (const auto& k : r.intervals)
				res.push_back(compute_estimate(
					r.funcName, r.func,
					r.funcExpected, k,
					r.tolerance, r.quiet,
					r.iterations
					));

			return res;
		}


		/// Register an equation for evaluation
		inline void equals(
			const std::string& name,
			Real evaluated, Real expected,
			Real tolerance = CHEBYSHEV_TOLERANCE, bool quiet = false) {

			equation_request r;
			r.funcName = name;
			r.evaluated = evaluated;
			r.expected = expected;
			r.tolerance = tolerance;
			r.quiet = quiet;

			state.equation_requests.push_back(r);
		}


		/// Register different equation evaluations
		inline void equals(
			const std::string& name,
			std::vector<std::array<Real, 2>> values,
			Real tolerance = CHEBYSHEV_TOLERANCE,
			bool quiet = false) {

			for (const auto& v : values)
				equals(name, v[0], v[1], tolerance, quiet);
		}


		/// Register a function for equation evaluation
		inline void equals(
			const std::string& name,
			RealFunction f,
			std::vector<std::array<Real, 2>> values,
			Real tolerance = CHEBYSHEV_TOLERANCE,
			bool quiet = false) {

			for (const auto& v : values)
				equals(name, f(v[0]), v[1], tolerance, quiet);
		}


		/// Test whether two real values are almost equal, to the given tolerance
		inline equation_result eval_equation(
			const std::string& name,
			Real evaluated, Real expected,
			Real tolerance = CHEBYSHEV_TOLERANCE, bool quiet = false) {

			equation_result eq;
			Real diff = std::abs(expected - evaluated);

			if(diff > tolerance) {
				state.failedTests++;
				eq.failed = true;
			} else {
			 	eq.failed = false;
			}

			state.totalTests++;

			eq.funcName = name;
			eq.diff = diff;
			eq.expected = expected;
			eq.evaluated = evaluated;
			eq.tolerance = tolerance;
			eq.quiet = quiet;

			state.equation_results[eq.funcName].push_back(eq);
			return eq;
		}


		inline equation_result eval_equation(equation_request er) {
			return eval_equation(er.funcName, er.evaluated, er.expected, er.tolerance, er.quiet);
		}


		/// Setup the precision testing environment
		inline void setup(std::string moduleName) {

			std::cout << "Starting precision testing of the " << moduleName << " module ..." << std::endl;
			state.moduleName = moduleName;
			state.failedTests = 0;
			state.totalTests = 0;

			if(state.outputToFile) {

				std::string filename;
				filename = state.outputFolder + "prec_" + moduleName + ".csv";

				if(state.outputFile.is_open())
					state.outputFile.close();

				state.outputFile.open(filename);

				if(!state.outputFile.is_open()) {
					std::cout << "Unable to open output file, results will NOT be saved!" << std::endl;
					state.outputToFile = false;
				}
			}
		}


		/// Run all requested error estimations and equation evaluations
		inline void run() {

			if(state.estimation_requests.size()) {

				if(!state.quiet) {
					std::cout << "\n" << std::left << std::setw(20) << "Function" << " | "
						<< std::setw(12) << "Int. Min." << " | "
						<< std::setw(12) << "Int. Max." << " | "
						<< std::setw(12) << "Mean Err." << " | "
						<< std::setw(12) << "RMS Err." << " | "
						<< std::setw(12)<< "Max Err." << " | "
						<< std::setw(12)<< "Rel. Err." << std::endl;
				}

				if(state.outputToFile)
					state.outputFile << "Function, Int. Min., Int. Max., Mean Err., "
						<< "RMS Err., Max Err., Rel. Err." << std::endl;

				for(const auto& r : state.estimation_requests) {
					
					auto res = compute_estimate(r);

					for(size_t i = 0; i < res.size(); i++) {

						if(state.estimateOnlyFailed && !res[i].failed)
							continue;

						if(!state.quiet) {

							std::cout << std::left << std::setw(20);

							if(i)	std::cout << "                    ";
							else	std::cout << res[i].funcName;
							
							std::cout << " | "
							<< std::setw(12) << res[i].k.a << " | "
							<< std::setw(12) << res[i].k.b << " | "
							<< std::setw(12) << res[i].mean_err << " | "
							<< std::setw(12) << res[i].rms_err << " | "
							<< std::setw(12) << res[i].max_err << " | "
							<< std::setw(12) << res[i].rel_err;

							if(res[i].failed)
								std::cout << "  FAILED";

							std::cout << std::endl;
						}

						if(state.outputToFile) {
							state.outputFile << res[i].funcName << ", " << res[i].k.a << ", " << res[i].k.b << ", "
								<< res[i].mean_err << ", " << res[i].rms_err << ", "
								<< res[i].max_err << ", " << res[i].rel_err << std::endl;
						}

					}

				}

				state.estimation_requests.clear();
			}


			if(state.equation_requests.size()) {

				if(!state.quiet) {
					std::cout << "\n" << std::setw(20) << "Function" << " | "
					 << std::setw(12) << "Eval. Value" << " | "
					 << std::setw(12) << "Exp. Value" << " | "
					 << std::setw(12) << "Diff." << " | "
					 << std::setw(12) << "Tol." << std::endl;
				}
				
				if(state.outputToFile)
					state.outputFile << "Function, Eval. Value, Exp. Value, Diff., Tol." << std::endl;

				for (size_t i = 0; i < state.equation_requests.size(); i++) {
					
					equation_result res = eval_equation(state.equation_requests[i]);

					if(state.equalsOnlyFailed && !res.failed)
						return;

					if(!state.quiet) {
						std::cout << std::setw(20);

						if(i) {
							if(state.equation_requests[i - 1].funcName == res.funcName)
								std::cout << "                    ";
							else
								std::cout << res.funcName;
						} else {
							std::cout << res.funcName;
						}

						std::cout << " | "
						<< std::setw(12) << res.evaluated << " | "
						<< std::setw(12) << res.expected << " | "
						<< std::setw(12) << res.diff << " | "
						<< std::setw(12) << res.tolerance;

						if(res.failed)
								std::cout << std::setw(8) << "  FAILED";

						std::cout << std::endl;
					}

					if(state.outputToFile) {
						state.outputFile << res.funcName << ", " << res.evaluated << ", " << res.expected
						<< ", " << res.diff << ", " << res.tolerance << std::endl;
					}
				}
			
				state.equation_requests.clear();
			}
		}


		/// Terminate precision testing
		inline void terminate(bool exit = true) {

			if(state.equation_requests.size() + state.equation_requests.size())
				run();

			std::cout << "\nFinished testing " << state.moduleName << std::endl;
			std::cout << state.totalTests << " total tests, " << state.failedTests << " failed (" <<
				(state.failedTests / (double) state.totalTests) * 100 << "%)" << std::endl;
				
			std::cout << "Results have been saved in "
				<< state.outputFolder << "prec_" << state.moduleName << ".csv" << std::endl;
			state.outputFile.close();
			
			state = prec_state();

			if(exit)
				std::exit(state.failedTests);
		}

	}

}
