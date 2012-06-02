# $ make
# $ make install (copy into bin directory)

CC = gcc
LINK = gcc
INSTALL = install

CFLAGS = -march=native -O3 -I../include -I.
LFLAGS = -march=native

GMPLIB = -L/usr/lib -lgmp 

BIN = bin
OBJ = obj
SRC = src

all: $(BIN)/newicktree $(BIN)/ranked_tree #$(BIN)/gmptest

# Generate object files
$(OBJ)/generate_sarray.o: $(SRC)/generate_sarray.c
	$(CC) $(CFLAGS) -o $@ -c $^

$(OBJ)/utils.o: $(SRC)/utils.c
	$(CC) $(CFLAGS) -o $@ -c $^

$(OBJ)/seqMain.o: $(SRC)/seqMain.c
	$(CC) $(CFLAGS) -o $@ -c $^

$(OBJ)/seqUtil.o: $(SRC)/seqUtil.c
	$(CC) $(CFLAGS) -o $@ -c $^

$(OBJ)/Newickform.o: $(SRC)/Newickform.c
	$(CC) $(CFLAGS) -o $@ -c $^

$(OBJ)/getopt.o: $(SRC)/getopt.c
	$(CC) $(CFLAGS) -o $@ -c $^

# Build binaries

#$(BIN)/gmptest:
#	gcc -L/usr/lib -lgmp gmp_test.c -o gmptest

$(BIN)/ranked_tree: $(OBJ)/Newickform.o $(OBJ)/generate_sarray.o $(OBJ)/utils.o $(OBJ)/getopt.o $(OBJ)/seqUtil.o
	$(LINK) -o $@ $^ $(LFLAGS)

$(BIN)/newicktree: $(OBJ)/seqMain.o $(OBJ)/seqUtil.o $(OBJ)/Newickform.o
	$(LINK) -o $@ $^ $(LFLAGS)

clean:
	-rm *.o $(BIN)/newicktree $(BIN)/ranked_tree
