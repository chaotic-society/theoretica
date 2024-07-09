default_target: all
.PHONY: all precision benchmark errors
all: precision benchmark errors

CXXFLAGS = -std=c++14 -I./src/ -Wall

precision:
	@echo Compiling \"precision\" example program ...
	@g++ examples/precision.cpp ${CXXFLAGS} -o ./precision

benchmark:
	@echo Compiling \"benchmark\" example program ...
	@g++ examples/benchmark.cpp ${CXXFLAGS} -o ./benchmark

errors:
	@echo Compiling \"errors\" example program ...
	@g++ examples/errors.cpp ${CXXFLAGS} -o ./errors
