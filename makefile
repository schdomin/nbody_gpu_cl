#ds compiler
CC = g++

#ds compiler flags
CFLAGS = -c -Wall

#ds default field
all: main

	$(CC) bin/CVector.o bin/CCubicDomain.o bin/main.o -o bin/nbody_cpu_cl

#ds object files
main:

	rm -rf bin
	mkdir bin
	$(CC) $(CFLAGS) src/CVector.cpp -o bin/CVector.o
	$(CC) $(CFLAGS) src/CCubicDomain.cpp -o bin/CCubicDomain.o
	$(CC) $(CFLAGS) src/main.cpp -o bin/main.o

#ds mark clean as independent
.PHONY: clean

#ds clean command
clean:

	rm -rf bin
