#ds compiler
CC = nvcc

#ds compiler flags
CFLAGS = -c -arch sm_13

#ds default field
all: main

	$(CC) bin/CVector.o bin/CCubicDomain.o bin/main.o -o bin/nbody_gpu_cl

#ds object files
main:

	rm -rf bin
	mkdir bin
	$(CC) $(CFLAGS) src/CVector.cu -o bin/CVector.o
	$(CC) $(CFLAGS) src/CCubicDomain.cu -o bin/CCubicDomain.o
	$(CC) $(CFLAGS) src/main.cu -o bin/main.o

#ds mark clean as independent
.PHONY: clean

#ds clean command
clean:

	rm -rf bin
