default_target: example

.PHONY: all example test

all: example test

example:
	g++ src/example.cpp -O2 -DUROBORO_INCLUDE_ALL -DUROBORO_X86 -o ./example

test:
	g++ test/test_real_analysis.cpp -DUROBORO_INCLUDE_ALL -DUROBORO_X86 -O2 -o test/test_real_analysis
	g++ test/test_algebra.cpp -DUROBORO_INCLUDE_ALL -DUROBORO_X86 -O2 -o test/test_algebra
	./test/test_real_analysis
	./test/test_algebra
