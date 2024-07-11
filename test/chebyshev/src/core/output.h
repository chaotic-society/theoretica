
///
/// @file output.h The output module, with formatting capabilities.
///

#ifndef CHEBYSHEV_OUTPUT_H
#define CHEBYSHEV_OUTPUT_H

#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <iomanip>

#include "../prec/prec_structures.h"
#include "../benchmark/benchmark_structures.h"
#include "../err/err_structures.h"


namespace chebyshev {

	/// @namespace output Functions to manage printing results.
	namespace output {

		/// @class field_options
		/// Custom options for printing a certain field.
		struct field_options {

			using FieldInterpreter = std::function<std::string(const std::string&)>;
			
			/// Width for the column associated with the field.
			unsigned int columnWidth = CHEBYSHEV_OUTPUT_WIDTH;

			/// A function which gets as input the value of a field as a string
			/// and returns a new string (e.g. "1" -> "FAIL" in the field "failed").
			FieldInterpreter fieldInterpreter = [](const std::string& s) { return s; };

			/// Additional custom options.
			std::map<std::string, long double> additionalFields {};

			// Default constructor
			field_options() {}

			// Construct field options from the custom column width.
			field_options(unsigned int columnWidth) : columnWidth(columnWidth) {}
		};


		/// @class output_state
		/// Global state of printing results to standard output.
		struct output_state {

			using OutputFormat_t = std::function<
			std::string(
				const std::vector<std::vector<std::string>>&,
				const std::vector<std::string>&,
				const output_state&)>;

			/// Relative or absolute path to output folder
			std::string outputFolder = "";

			/// Map of field name to output string
			/// (e.g. "maxErr" -> "Max Err.").
			std::map<std::string, std::string> fieldNames {};

			/// Options for the different fields.
			std::map<std::string, field_options> fieldOptions {};

			/// The output files, indexed by filename
			std::map<std::string, std::ofstream> outputFiles {};

			/// Default width for a field.
			unsigned int defaultColumnWidth = CHEBYSHEV_OUTPUT_WIDTH;

			/// The number of digits to show in scientific notation.
			unsigned int outputPrecision = 1;

			/// The output format to use to print to standard output.
			OutputFormat_t outputFormat {};

			/// The default output format to use for files,
			/// when no format has been set for a file.
			OutputFormat_t defaultFileOutputFormat {};

			/// The output format to use for a specific file, by filename.
			std::map<std::string, OutputFormat_t> fileOutputFormat {};

			/// Whether to output to standard output.
			bool quiet = false;

			/// Whether the output module was setup.
			bool wasSetup = false;

		} state;


		/// A function which converts the table entries of a row
		/// to a string to print (e.g. adding separators and padding).
		/// @see output_state::OutputFormat_t
		using OutputFormat = output_state::OutputFormat_t;


		/// @namespace format Output formatting functions.
		namespace format {
			

			/// Bare bone output format which just prints the result
			/// table as is, without any formatting beyond adjusting column width.
			inline OutputFormat barebone() {

				return [](
					const std::vector<std::vector<std::string>>& table,
					const std::vector<std::string>& fields,
					const output_state& state) -> std::string {

					if(!table.size())
						return "";

					std::stringstream result;

					for (size_t i = 0; i < table.size(); ++i) {

						if(table[i].size() != fields.size()) {
							std::runtime_error(
								"Number of columns and fields argument must have "
								"the same size in output::format::barebone");
						}

						for (size_t j = 0; j < table[i].size(); ++j) {

							auto it = state.fieldOptions.find(fields[j]);

							if(it != state.fieldOptions.end() && i)
								result << std::setw(it->second.columnWidth)
								<< std::left << it->second.fieldInterpreter(table[i][j]);
							else if(it != state.fieldOptions.end())
								result << std::setw(it->second.columnWidth)
								<< std::left << table[i][j];
							else
								result << std::setw(state.defaultColumnWidth)
								<< std::left << table[i][j];
						}

						result << "\n";
					}

					return result.str();
				};
			}


			/// Simple output format which prints the fields
			/// separated by the separator string and padding, if enabled.
			/// The OutputFormat is returned as a lambda function.
			/// This format is a good starting point if you want to implement
			/// your own custom output format.
			inline OutputFormat simple() {

				return [](
					const std::vector<std::vector<std::string>>& table,
					const std::vector<std::string>& fields,
					const output_state& state) -> std::string {

					if(!table.size())
						return "";

					std::stringstream result;
					std::stringstream header_str;

					header_str << " | ";

					for (size_t i = 0; i < table[0].size(); ++i) {

						auto it = state.fieldOptions.find(fields[i]);

						if(it != state.fieldOptions.end())
							header_str << std::setw(it->second.columnWidth) << table[0][i] << " | ";
						else
							header_str << std::setw(state.defaultColumnWidth) << table[0][i] << " | ";
					}

					std::string header = header_str.str();
					std::string decoration = " +";
					for (size_t i = 4; i < header.size(); ++i)
						decoration += "-";
					decoration += "+ \n";

					for (size_t i = 1; i < table.size(); ++i) {

						if(table[i].size() != fields.size()) {
							std::runtime_error(
								"Number of columns and <fields> argument must have "
								"the same size in output::format::simple");
						}

						result << " | ";

						for (size_t j = 0; j < table[i].size(); ++j) {

							auto it = state.fieldOptions.find(fields[j]);

							if(it != state.fieldOptions.end())
								result << std::setw(it->second.columnWidth)
								<< it->second.fieldInterpreter(table[i][j]) << " | ";
							else
								result << std::setw(state.defaultColumnWidth)
								<< table[i][j] << " | ";
						}

						result << "\n";
					}

					return decoration
						+ header + "\n"
						+ decoration
						+ result.str()
						+ decoration;
				};
			}


			/// Fancy output format which uses Unicode characters
			/// to print a continuous outline around the table.
			/// The OutputFormat is returned as a lambda function.
			inline OutputFormat fancy() {

				return [](
					const std::vector<std::vector<std::string>>& table,
					const std::vector<std::string>& fields,
					const output_state& state) -> std::string {

					if(!table.size())
						return "";

					// Effective length of the string
					// (needed because Unicode is used)
					size_t eff_length = 0;
					std::stringstream header_str;
					
					header_str << " │ ";
					eff_length += 3;

					for (size_t i = 0; i < table[0].size(); ++i) {

						auto it = state.fieldOptions.find(fields[i]);

						if(it != state.fieldOptions.end()) {
							header_str << std::setw(it->second.columnWidth)
							<< table[0][i] << " │ ";
							eff_length += it->second.columnWidth;
						} else {
							header_str << std::setw(state.defaultColumnWidth)
							<< table[0][i] << " │ ";
							eff_length += state.defaultColumnWidth;
						}

						eff_length += 3;
					}

					std::string header = " ┌";

					// Upper outline
					for (size_t i = 4; i < eff_length; ++i)
						header += "─";
					header += "┐ \n";

					header += header_str.str() + "\n";

					// Lower outline
					header += " ├";
					for (size_t i = 4; i < eff_length; ++i)
						header += "─";
					header += "┤ \n";

					std::stringstream result;

					for (size_t i = 1; i < table.size(); ++i) {

						if(table[i].size() != fields.size()) {
							std::runtime_error(
								"Number of columns and <fields> argument must have "
								"the same size in output::format::fancy");
						}

						result << " │ ";

						for (size_t j = 0; j < table[i].size(); ++j) {

							auto it = state.fieldOptions.find(fields[j]);

							if(it != state.fieldOptions.end())
								result << std::setw(it->second.columnWidth)
								<< it->second.fieldInterpreter(table[i][j]) << " │ ";
							else
								result << std::setw(state.defaultColumnWidth)
								<< table[i][j] << " │ ";
						}

						result << "\n";
					}

					std::string underline = " └";
					for (size_t i = 4; i < eff_length; ++i) {
						underline += "─";
					}
					underline += "┘ \n";

					return header + result.str() + underline;
				};
			}


			/// Format function for CSV format files.
			/// The OutputFormat is returned as a lambda function.
			///
			/// @param separator The string to print between
			/// different fields (defaults to ",").
			inline OutputFormat csv(const std::string& separator = ",") {

				return [separator](
					const std::vector<std::vector<std::string>>& table,
					const std::vector<std::string>& fields,
					const output_state& state) -> std::string {

					std::stringstream s;

					for (size_t i = 0; i < table.size(); ++i) {

						if(table[i].size() != fields.size()) {
							throw std::runtime_error(
								"Number of columns and <fields> argument must have "
								"the same size in output::format::csv");
						}


						for (size_t j = 0; j < table[i].size(); ++j) {

							auto it = state.fieldOptions.find(fields[j]);

							if(it != state.fieldOptions.end() && i)
								s << "\"" << it->second.fieldInterpreter(table[i][j]) << "\"";
							else
								s << "\"" << table[i][j] << "\"";

							if(j != table[i].size() - 1)
								s << separator;
						}

						s << "\n";
					}

					return s.str();
				};
			}


			/// Format the table as Markdown.
			/// The OutputFormat is returned as a lambda function.
			inline OutputFormat markdown() {

				return [](
					const std::vector<std::vector<std::string>>& table,
					const std::vector<std::string>& fields,
					const output_state& state) -> std::string {

					if(!table.size())
						return "";

					std::stringstream result;
					std::stringstream header_str;

					header_str << "|";

					for (size_t i = 0; i < table[0].size(); ++i) {

						auto it = state.fieldOptions.find(fields[i]);

						if(it != state.fieldOptions.end())
							header_str << std::setw(it->second.columnWidth)
							<< std::left << table[0][i];
						else
							header_str << std::setw(state.defaultColumnWidth)
							<< std::left << table[0][i];

						header_str << "|";
					}

					std::string header = header_str.str();
					std::string decoration = "|";
					for (size_t i = 1; i < header.size() - 1; ++i)
						decoration += (header[i] == '|') ? "|" : "-";
					decoration += "|\n";

					for (size_t i = 1; i < table.size(); ++i) {

						if(table[i].size() != fields.size()) {
							std::runtime_error(
								"Number of columns and <fields> argument must have "
								"the same size in output::format::markdown");
						}

						result << "|";

						for (size_t j = 0; j < table[i].size(); ++j) {

							auto it = state.fieldOptions.find(fields[j]);

							if(it != state.fieldOptions.end())
								result << std::setw(it->second.columnWidth)
								<< std::left << it->second.fieldInterpreter(table[i][j]);
							else
								result << std::setw(state.defaultColumnWidth)
								<< std::left << table[i][j];

							result << "|";
						}

						result << "\n";
					}

					return header + "\n" + decoration + result.str();
				};
			}


			/// Format the table as a LaTeX table in the tabular environment.
			/// The OutputFormat is returned as a lambda function.
			inline OutputFormat latex() {

				return [=](
					const std::vector<std::vector<std::string>>& table,
					const std::vector<std::string>& fields,
					const output_state& state) -> std::string {

					if(!table.size())
						return "";

					std::stringstream result;
					result << "\\begin{tabular}{";

					if(fields.size())
						result << "|";

					for (unsigned int i = 0; i < fields.size(); ++i)
						result << "c|";

					result << "}\n\\hline\n";

					for (size_t i = 0; i < table[0].size(); ++i) {

						result << table[0][i];

						if(i != table[0].size() - 1)
							result << " & ";
					}
					result << " \\\\\n\\hline\n";

					for (size_t i = 1; i < table.size(); ++i) {

						if(table[i].size() != fields.size()) {
							std::runtime_error(
								"Number of columns and <fields> argument must have "
								"the same size in output::format::latex");
						}

						for (size_t j = 0; j < table[i].size(); ++j) {

							auto it = state.fieldOptions.find(fields[j]);

							if(it != state.fieldOptions.end())
								result << it->second.fieldInterpreter(table[i][j]);
							else
								result << table[i][j];

							if(j != table[i].size() - 1)
								result << " & ";
						}

						result << " \\\\\n";
					}

					result << "\\hline\n\\end{tabular}\n";

					return result.str();
				};
			}

		}


		/// Setup printing to the output stream with default options.
		inline void setup() {

			// Estimate fields
			state.fieldNames["funcName"] = "Function";
			state.fieldNames["maxErr"] = "Max Err.";
			state.fieldNames["meanErr"] = "Mean Err.";
			state.fieldNames["rmsErr"] = "RMS Err.";
			state.fieldNames["relErr"] = "Rel. Err.";
			state.fieldNames["absErr"] = "Abs. Err.";
			state.fieldNames["tolerance"] = "Tolerance";
			state.fieldNames["failed"] = "Result";
			state.fieldNames["iterations"] = "Iterations";

			// Equation fields
			state.fieldNames["difference"] = "Difference";
			state.fieldNames["evaluated"] = "Evaluated";
			state.fieldNames["expected"] = "Expected";

			// Benchmark fields
			state.fieldNames["totalRuntime"] = "Tot. Time (ms)";
			state.fieldNames["averageRuntime"] = "Avg. Time (ms)";
			state.fieldNames["runsPerSecond"] = "Runs per Sec.";
			state.fieldNames["runs"] = "Runs";

			// Error checking
			state.fieldNames["correctType"] = "Correct Type";
			state.fieldNames["description"] = "Description";
			state.fieldNames["expectedFlags"] = "Exp. Flags";
			state.fieldNames["thrown"] = "Has Thrown";

			// Set wider column width for some fields
			state.fieldOptions["funcName"].columnWidth = 16;
			state.fieldOptions["averageRuntime"].columnWidth = 14;
			state.fieldOptions["runsPerSecond"].columnWidth = 14;
			state.fieldOptions["description"].columnWidth = 20;

			// Set a special field interpreter for the "failed" field
			state.fieldOptions["failed"].fieldInterpreter = [](const std::string& s) {
				if(s == "0") return "PASS";
				else if(s == "1") return "FAIL";
				else return "UNKNOWN";
			};

			// Set the default output formats
			state.outputFormat = format::fancy();
			state.defaultFileOutputFormat = format::csv();

			state.wasSetup = true;
		}


		/// Terminate the output module by closing all output files
		/// and resetting its state.
		inline void terminate() {

			for (auto& p : state.outputFiles)
				if(p.second.is_open())
					p.second.close();

			state.outputFiles.clear();
			state = output_state();
		}


		/// Resolve the field of an estimate result by name,
		/// returning the value as a string.
		///
		/// @param fieldName The name of the field to resolve
		/// @param r The estimate result to read the fields of
		inline std::string resolve_field(
			const std::string& fieldName, prec::estimate_result r) {

			std::stringstream value;

			if(fieldName == "funcName") {
				value << r.funcName;
			} else if(fieldName == "maxErr") {
				value << std::scientific
					<< std::setprecision(state.outputPrecision)
					<< r.maxErr;
			} else if(fieldName == "meanErr") {
				value << std::scientific
					<< std::setprecision(state.outputPrecision)
					<< r.meanErr;
			} else if(fieldName == "rmsErr") {
				value << std::scientific
					<< std::setprecision(state.outputPrecision)
					<< r.rmsErr;
			} else if(fieldName == "relErr") {
				value << std::scientific
					<< std::setprecision(state.outputPrecision)
					<< r.relErr;
			} else if(fieldName == "absErr") {
				value << std::scientific
					<< std::setprecision(state.outputPrecision)
					<< r.absErr;
			} else if(fieldName == "tolerance") {
				value << std::scientific
					<< std::setprecision(state.outputPrecision)
					<< r.tolerance;
			} else if(fieldName == "failed") {
				value << r.failed;
			} else {
				
				if(r.additionalFields.find(fieldName) == r.additionalFields.end())
					return "";

				value << r.additionalFields[fieldName];
			}

			return value.str();
		}


		/// Resolve the field of an equation result by name,
		/// returning the value as a string.
		///
		/// @param fieldName The name of the field to resolve
		/// @param r The equation result to read the fields of
		inline std::string resolve_field(
			const std::string& fieldName, prec::equation_result r) {

			std::stringstream value;

			if(fieldName == "funcName") {
				value << r.funcName;
			} else if(fieldName == "evaluated") {
				value << r.evaluated;
			} else if(fieldName == "expected") {
				value << r.expected;
			} else if(fieldName == "difference") {
				value << std::scientific
					<< std::setprecision(state.outputPrecision)
					<< r.difference;
			} else if(fieldName == "tolerance") {
				value << std::scientific
					<< std::setprecision(state.outputPrecision)
					<< r.tolerance;
			} else if(fieldName == "failed") {
				value << r.failed;
			} else {
				
				if(r.additionalFields.find(fieldName) == r.additionalFields.end())
					return "";

				value << r.additionalFields[fieldName];
			}

			return value.str();
		}


		/// Resolve the field of a benchmark result by name,
		/// returning the value as a string.
		inline std::string resolve_field(
			const std::string& fieldName, benchmark::benchmark_result r) {

			std::stringstream value;

			if(fieldName == "funcName") {
				value << r.funcName;
			} else if(fieldName == "runs") {
				value << r.runs;
			} else if(fieldName == "iterations") {
				value << r.iterations;
			} else if(fieldName == "totalRuntime") {
				value << r.totalRuntime;
			} else if(fieldName == "averageRuntime") {
				value << r.averageRuntime;
			} else if(fieldName == "runsPerSecond") {
				if(r.runsPerSecond > 1000)
					value << uint64_t(r.runsPerSecond);
				else
					value << r.runsPerSecond;
			} else if(fieldName == "failed") {
				value << r.failed;
			} else {
				
				if(r.additionalFields.find(fieldName) == r.additionalFields.end())
					return "";

				value << r.additionalFields[fieldName];
			}

			return value.str();
		}


		/// Resolve the field of an assertion result by name,
		/// returning the value as a string.
		///
		/// @param fieldName The name of the field to resolve
		/// @param r The assertion result to read the fields of
		inline std::string resolve_field(
			const std::string& fieldName, err::assert_result r) {

			std::stringstream value;

			if(fieldName == "funcName") {
				value << r.funcName;
			} else if(fieldName == "evaluated") {
				value << r.evaluated;
			} else if(fieldName == "description") {
				value << r.description;
			} else if(fieldName == "failed") {
				value << r.failed;
			} else {
				return "";
			}

			return value.str();
		}


		/// Resolve the field of an errno checking result by name,
		/// returning the value as a string.
		///
		/// @param fieldName The name of the field to resolve
		/// @param r The errno checking result to read the fields of
		inline std::string resolve_field(
			const std::string& fieldName, err::errno_result r) {

			std::stringstream value;

			if(fieldName == "funcName") {
				value << r.funcName;
			} else if(fieldName == "evaluated") {
				value << r.evaluated;
			} else if(fieldName == "expectedFlags") {

				int res_flag = 0xFFFFFFFF;
				for (int flag : r.expectedFlags)
					res_flag &= flag;

				value << res_flag;

			} else if(fieldName == "failed") {
				value << r.failed;
			} else {
				return "";
			}

			return value.str();
		}


		/// Resolve the field of an exception checking result by name,
		/// returning the value as a string.
		///
		/// @param fieldName The name of the field to resolve
		/// @param r The exception checking result to read the fields of
		inline std::string resolve_field(
			const std::string& fieldName, err::exception_result r) {

			std::stringstream value;

			if(fieldName == "funcName") {
				value << r.funcName;
			} else if(fieldName == "thrown") {
				value << r.thrown;
			} else if(fieldName == "correctType") {
				value << r.correctType;
			} else if(fieldName == "failed") {
				value << r.failed;
			} else {
				return "";
			}

			return value.str();
		}


		/// Generate a table of results as a string matrix to pass to
		/// a specific formatter of OutputFormat type.
		///
		/// @param results The map of test results of any type
		/// @param fields The fields of the test results to write
		/// to each column, in order.
		/// @return A string matrix representing the results as a table.
		template<typename ResultType>
		inline auto generate_table(
			const std::map<std::string, std::vector<ResultType>>& results,
			const std::vector<std::string>& fields) {

			std::vector<std::vector<std::string>> table;

			// Construct header
			std::vector<std::string> header (fields.size());
			for (size_t i = 0; i < fields.size(); ++i) {

				const auto it = state.fieldNames.find(fields[i]);

				// Associate string to field name
				if(it != state.fieldNames.end())
					header[i] = it->second;
				else
					header[i] = fields[i];
			}
			table.emplace_back(header);

			// Construct rows
			for (const auto& p : results) {
				for (const auto& result : p.second) {

					std::vector<std::string> row (fields.size());

					for (size_t i = 0; i < fields.size(); ++i)
						row[i] = resolve_field(fields[i], result);

					table.emplace_back(row);
				}
			}

			return table;
		}


		/// Print the test results to standard output and output files
		/// with their given formats, defaulting to state.outputFiles
		/// if no filenames are specified.
		///
		/// @param results The map of test results, of any type.
		/// @param fields The fields of the test results to write, in order.
		/// @param filenames The names of the output files
		/// (if this list is empty, state.outputFiles will be used).
		template<typename ResultType>
		inline void print_results(
			const std::map<std::string, std::vector<ResultType>>& results,
			const std::vector<std::string>& fields,
			const std::vector<std::string>& filenames) {

			// Skip output on no test case results
			if(results.empty())
				return;

			// Table data as a string matrix
			std::vector<std::vector<std::string>> table = generate_table(results, fields);

			// Write to standard output
			if(!state.quiet)
				std::cout << "\n" << state.outputFormat(table, fields, state) << "\n";

			// Write to the specified filenames or to the
			// default output files if no filename is specified.
			if(filenames.size()) {

				for (size_t i = 0; i < filenames.size(); ++i) {

					std::ofstream file (filenames[i]);
					if(!file.is_open()) {
						std::cout << "Cannot write to output file: " << filenames[i] << std::endl;
						continue;
					}

					// Apply formatting according to set options
					auto it = state.fileOutputFormat.find(filenames[i]);

					if(it != state.fileOutputFormat.end())
						file << it->second(table, fields, state);
					else
						file << state.defaultFileOutputFormat(table, fields, state);

					std::cout << "Results have been saved in: " << filenames[i] << std::endl;
				}

			} else {

				for (auto& p : state.outputFiles) {

					if(!p.second.is_open()) {
						std::cout << "Cannot write to output file: " << p.first << std::endl;
						continue;
					}
					
					// Apply formatting according to set options
					auto it = state.fileOutputFormat.find(p.first);

					if(it != state.fileOutputFormat.end())
						p.second << it->second(table, fields, state);
					else
						p.second << state.defaultFileOutputFormat(table, fields, state);

					std::cout << "Results have been saved in: " << p.first << std::endl;
				}

			}

		}
	}
}

#endif
