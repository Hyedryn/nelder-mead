
FLAGS = -std=c99 -Wall -Wextra -pedantic
CC_OPT = -O3 -ffast-math -fno-common -funroll-loops
BIN = nm
SRC = main.c nelder_mead.c point.c mfobj.c

compile:
	gcc $(FLAGS) $(CC_OPT) -o $(BIN) $(SRC) -lm

test: compile
	time ./$(BIN) -2.10 -3.04 4.50 
