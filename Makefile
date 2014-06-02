CC=mpicc
CFLAGS=-Wall -std=gnu99 -DDEBUG -ggdb
sor: main.o test.o 
	$(CC) $(CFLAGS) -o sor main.o test.o
