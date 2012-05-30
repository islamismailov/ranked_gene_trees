# $ make

CC = gcc
LINK = gcc
INSTALL = install
CFLAGS = -march=native -O3 -I../include -I.
LFLAGS = -march=native 

all: gmptest newicktree spec_array.o utils.o

spec_array.o: generate_sarray.c
	gcc -c generate_sarray.c

utils.o: utils.c
	gcc -c utils.c

gmptest:
	gcc -L/usr/lib -lgmp gmp_test.c -o gmptest

newicktree: seqMain.o seqUtil.o Newickform.o
	gcc -o newicktree seqMain.o seqUtil.o Newickform.o

seqMain.o: seqMain.c
	gcc -c seqMain.c

seqUtil.o: seqUtil.c
	gcc -c seqUtil.c

Newickform.o: Newickform.c
	gcc -c Newickform.c

clean:
	-rm *.o newicktree gmptest
