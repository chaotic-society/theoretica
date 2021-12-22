default_target: example

.PHONY: all example test

all: example test

example:
	g++ src/example.cpp -O2 -DUROBORO_INCLUDE_ALL -DUROBORO_X86 -o ./example

test:
	g++ test/test.cpp -DUROBORO_INCLUDE_ALL -DUROBORO_X86 -O2 -o test/test
	./test/test
