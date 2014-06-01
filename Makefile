CC=mpicc
CFLAGS=-Wall -std=gnu99 -DDEBUG -ggdb
sor: main.c
	$(CC) $(CFLAGS) -o sor main.c
