# CaRaSh local spherical shell experiment (GNU Make)

CC = gcc
CFLAGS = -std=c99 -O2 -Wall -Wextra
LDFLAGS = -lm

# Default lattice side (odd integer)
L ?= 17

all: carash

carash: carash.c
	$(CC) $(CFLAGS) -DL=$(L) -o $@ $< $(LDFLAGS)

run: carash
	./carash

clean:
	rm -f carash carash31 *.o

.PHONY: all run clean
