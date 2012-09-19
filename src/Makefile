CC = gcc -g
CFLAGS = -O2 -Wall

sega: main.o train.o seq2qu.o util.o
	$(CC) $(CFLAGS) -o ../bin/sega main.o train.o seq2qu.o util.o -lm -lz
main.o: main.c
	$(CC) $(CFLAGS) -c main.c
train.o: train.c util.h
	$(CC) $(CFLAGS) -c train.c
seq2qu.o: seq2qu.c util.h
	$(CC) $(CFLAGS) -c seq2qu.c
util.o: util.c util.h
	$(CC) $(CFLAGS) -c util.c
