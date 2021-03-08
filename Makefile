default_target: all

all:
	g++ src/main.cpp -O2 -o ./main

labcalc:
	g++ tools/labcalc.cpp -O2 -I./src/ -o ./labcalc
