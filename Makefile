default_target: example

.PHONY: all example test labcalc

all: example test labcalc

example:
	g++ src/example.cpp -O2 -o ./example

test:
	g++ test/test.cpp -o test/test
	./test/test
