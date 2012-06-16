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

all: ensure_dirs $(BIN)/newicktree $(BIN)/ranked_tree

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
$(BIN)/ranked_tree: $(OBJ)/Newickform.o $(OBJ)/generate_sarray.o $(OBJ)/utils.o $(OBJ)/getopt.o $(OBJ)/seqUtil.o
	$(LINK) $(LFLAGS) -o $@ $^

$(BIN)/newicktree: $(OBJ)/seqMain.o $(OBJ)/seqUtil.o $(OBJ)/Newickform.o
	$(LINK) $(LFLAGS) -o $@ $^

.PHONY: ensure_dirs clean

ensure_dirs:
	@if [ ! -d $(BIN) ]; then mkdir $(BIN); fi
	@if [ ! -d $(OBJ) ]; then mkdir $(OBJ); fi

clean:
	rm -rf $(OBJ)/*.o $(BIN)/newicktree $(BIN)/ranked_tree
