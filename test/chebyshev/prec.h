
///
/// @file prec.h Precision estimation of real functions
///

#pragma once

#include "core/common.h"
#include "core/interval.h"
#include "core/prec_def.h"

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

			/// Default fail function
			FailFunction defaultFailFunction = fail_on_max_err;

			/// Default tolerance on max absolute error
			Real defaultTolerance = CHEBYSHEV_TOLERANCE;

			/// Recorded estimation requests
			std::vector<estimate_request> estimationRequests;

			/// Recorded custom estimation requests
			std::vector<estimate_custom_request> estimationCustomRequests;

			/// Recorded equation requests
			std::vector<equation_request> equationRequests;

			/// Results of precision testing
			std::map<std::string, std::vector<estimate_result>> estimationResults;
			
			/// Results of equations
			std::map<std::string, std::vector<equation_result>> equationResults;

			/// Target tests marked for execution
			/// (all tests will be executed if empty)
			std::map<std::string, bool> pickedTests;

		} state;


		/// Register a function for error estimation
		inline void estimate(
			const std::string& name,
			RealFunction fApprox,
			RealFunction fExp,
			interval k,
			Real tolerance = state.defaultTolerance,
			bool quiet = false,
			unsigned int n = state.defaultIterations,
			FailFunction fail = state.defaultFailFunction) {

			estimate_request r;
			r.funcName = name;
			r.func = fApprox;
			r.funcExpected = fExp;
			r.intervals.push_back(k);
			r.tolerance = tolerance;
			r.quiet = quiet;
			r.iterations = n;
			r.fail = fail;

			state.estimationRequests.push_back(r);
		}


		/// Register a function for error estimation on multiple intervals
		inline void estimate(
			const std::string& name,
			RealFunction fApprox,
			RealFunction fExp,
			std::vector<interval> intervals,
			Real tolerance = state.defaultTolerance,
			bool quiet = false,
			unsigned int n = state.defaultIterations,
			FailFunction fail = state.defaultFailFunction) {

			estimate_request r;
			r.funcName = name;
			r.func = fApprox;
			r.funcExpected = fExp;
			r.intervals = intervals;
			r.tolerance = tolerance;
			r.quiet = quiet;
			r.iterations = n;
			r.fail = fail;

			state.estimationRequests.push_back(r);
		}


		inline void estimate(
			const std::string& name,
			CustomEstimateFunction f,
			interval k,
			Real tolerance = state.defaultTolerance,
			bool quiet = false,
			unsigned int n = state.defaultIterations,
			FailFunction fail = state.defaultFailFunction) {

			estimate_custom_request r;
			r.funcName = name;
			r.f = f;
			r.intervals.push_back(k);
			r.tolerance = tolerance;
			r.quiet = quiet;
			r.iterations = n;

			state.estimationCustomRequests.push_back(r);
		}


		inline void estimate(
			const std::string& name,
			CustomEstimateFunction f,
			std::vector<interval> intervals,
			Real tolerance = state.defaultTolerance,
			bool quiet = false,
			unsigned int n = state.defaultIterations,
			FailFunction fail = state.defaultFailFunction) {

			estimate_custom_request r;
			r.funcName = name;
			r.f = f;
			r.intervals = intervals;
			r.tolerance = tolerance;
			r.quiet = quiet;
			r.iterations = n;

			state.estimationCustomRequests.push_back(r);
		}


		/// Estimate the precision of a real function on a single interval
		inline estimate_result compute_estimate(
			std::string name,
			RealFunction fApprox,
			RealFunction fExp,
			interval k,
			Real tolerance = state.defaultTolerance,
			bool quiet = false,
			unsigned int n = state.defaultIterations,
			FailFunction fail = state.defaultFailFunction) {

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
			result.rel_err = std::abs((sum * dx / 3.0) / (sum_abs * dx / 3.0));
			result.tolerance = tolerance;

			result.failed = fail(result);

			state.estimationResults[result.funcName].push_back(result);
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
					r.iterations, r.fail
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

			state.equationRequests.push_back(r);
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

			state.equationResults[eq.funcName].push_back(eq);
			return eq;
		}


		inline equation_result eval_equation(equation_request er) {
			return eval_equation(er.funcName, er.evaluated, er.expected, er.tolerance, er.quiet);
		}


		/// Setup the precision testing environment
		inline void setup(std::string moduleName, int argc = 0, const char** argv = nullptr) {


			// Initialize pick list
			if(argc && argv) {
				for (int i = 1; i < argc; ++i) {
					state.pickedTests[argv[i]] = true;
				}
			}

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


		/// Print an estimation result
		inline void print_estimate(const estimate_result& r, size_t i = 0) {

			std::cout << std::left << std::setw(20);

			if(i)	std::cout << "                    ";
			else	std::cout << r.funcName;
			
			std::cout << " | "
			<< std::setw(12) << r.k.a << " | "
			<< std::setw(12) << r.k.b << " | "
			<< std::setw(12) << r.mean_err << " | "
			<< std::setw(12) << r.rms_err << " | "
			<< std::setw(12) << r.max_err << " | "
			<< std::setw(12) << r.rel_err;

			if(r.failed)
				std::cout << "  FAILED";

			std::cout << std::endl;
		}


		/// Print an equation result
		inline void print_equation(const equation_result& r, size_t i = 0) {

			std::cout << std::setw(20);

			if(i) {
				if(state.equationRequests[i - 1].funcName == r.funcName)
					std::cout << "                    ";
				else
					std::cout << r.funcName;
			} else {
				std::cout << r.funcName;
			}

			std::cout << " | "
			<< std::setw(12) << r.evaluated << " | "
			<< std::setw(12) << r.expected << " | "
			<< std::setw(12) << r.diff << " | "
			<< std::setw(12) << r.tolerance;

			if(r.failed)
					std::cout << std::setw(8) << "  FAILED";

			std::cout << std::endl;
		}


		/// Run all requested error estimations and equation evaluations
		inline void run() {

			if(state.estimationRequests.size() + state.estimationCustomRequests.size()) {

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

				for(const auto& r : state.estimationRequests) {

					if(!state.pickedTests.empty() && !state.pickedTests[r.funcName])
						continue;
					
					auto res = compute_estimate(r);

					for(size_t i = 0; i < res.size(); i++) {

						if(res[i].failed)
							state.failedTests++;

						// Skip test if only picked tests are to be executed
						if(state.estimateOnlyFailed && !res[i].failed)
							continue;

						if(!state.quiet)
							print_estimate(res[i], i);

						if(state.outputToFile) {
							state.outputFile << res[i].funcName << ", " << res[i].k.a << ", " << res[i].k.b << ", "
								<< res[i].mean_err << ", " << res[i].rms_err << ", "
								<< res[i].max_err << ", " << res[i].rel_err << std::endl;
						}

					}

				}
				state.estimationRequests.clear();

				for(const auto& r : state.estimationCustomRequests) {

					// Skip test if only picked tests are to be executed
					if(!state.pickedTests.empty() && !state.pickedTests[r.funcName])
						continue;
					
					std::vector<estimate_result> res;

					for (size_t j = 0; j < r.intervals.size(); ++j) {

						res.push_back(r.f(r.intervals[j], r.tolerance, r.iterations));
						res[j].funcName = r.funcName;
						res[j].k = r.intervals[j];
						res[j].tolerance = r.tolerance;
						res[j].iterations = r.iterations;
						res[j].quiet = r.quiet;

						state.totalTests++;

						if(res[j].failed)
							state.failedTests++;
					}

					for(size_t i = 0; i < res.size(); i++) {

						if(state.estimateOnlyFailed && !res[i].failed)
							continue;

						if(!state.quiet)
							print_estimate(res[i], i);

						if(state.outputToFile) {
							state.outputFile << res[i].funcName << ", " << res[i].k.a << ", " << res[i].k.b << ", "
								<< res[i].mean_err << ", " << res[i].rms_err << ", "
								<< res[i].max_err << ", " << res[i].rel_err << std::endl;
						}

					}

				}
				state.estimationCustomRequests.clear();
			}


			if(state.equationRequests.size()) {

				if(!state.quiet) {
					std::cout << "\n" << std::setw(20) << "Function" << " | "
					 << std::setw(12) << "Eval. Value" << " | "
					 << std::setw(12) << "Exp. Value" << " | "
					 << std::setw(12) << "Diff." << " | "
					 << std::setw(12) << "Tol." << std::endl;
				}
				
				if(state.outputToFile)
					state.outputFile << "Function, Eval. Value, Exp. Value, Diff., Tol." << std::endl;

				for (size_t i = 0; i < state.equationRequests.size(); i++) {

					// Skip test if only picked tests are to be executed
					if( !state.pickedTests.empty() &&
						!state.pickedTests[state.equationRequests[i].funcName])
						continue;
					
					equation_result res = eval_equation(state.equationRequests[i]);

					if(state.equalsOnlyFailed && !res.failed)
						return;

					if(!state.quiet)
						print_equation(res, i);

					if(state.outputToFile) {
						state.outputFile << res.funcName << ", " << res.evaluated << ", " << res.expected
						<< ", " << res.diff << ", " << res.tolerance << std::endl;
					}
				}
			
				state.equationRequests.clear();
			}
		}


		/// Terminate precision testing
		inline void terminate(bool exit = true) {

			if(state.estimationRequests.size() +
				state.equationRequests.size() +
				state.estimationCustomRequests.size())
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
