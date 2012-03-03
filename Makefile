# PREDEFINED VARIABLES
CFLAGS = -O3 -Wall -I. -Wno-unused
CC = gcc

ifeq ($(CC),gcc)
	LFLAGS = -lm
endif


#GLOBAL VARIABLES
PARAMFILE = param.txt
OUTPUTFILE = output

# MODULES
MODULES = io_lib initialize observables rk4 rk45 sources
MODULESC = $(MODULES:=.c)
MODULESH = $(MODULES:=.h)
MODULESO = $(MODULES:=.o)
PROGRAM = main
EXECUTABLE = exec

#Use python
first: clean $(EXECUTABLE).out
	python write_param.py

# LINKING
$(EXECUTABLE).out: $(PROGRAM).o $(MODULESO)
	$(CC) -o $@ $^ $(LFLAGS)
	rm -rf $(MODULESO) $(PROGRAM).o

# MODULE'S RULE
$(MODULES): $(MODULESC) $(MODULESH)
	$(CC) $@.c $(CFLAGS) -o $@.o

# PROGRAM COMPILATION
$(PROGRAM).o: $(PROGRAM).c $(MODULESH)
	$(CC) $(CFLAGS)   -c -o $@ $<

run:
	./$(EXECUTABLE).out $(PARAMFILE) $(OUTPUTFILE) 

clean:
	rm -rf *.out *.o *.c~ *.h~

all: clean $(EXECUTABLE).out
	./$(EXECUTABLE).out $(PARAMFILE) $(OUTPUTFILE) 


