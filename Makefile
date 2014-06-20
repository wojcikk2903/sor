CC=mpicc
CFLAGS=-Wall -std=gnu99 -ggdb
sor: main.o test.o io.o
	$(CC) $(CFLAGS) -o sor main.o test.o io.o
