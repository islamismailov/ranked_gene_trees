# $ make
# $ make install (copy into bin directory)

CC = gcc
LINK = gcc
INSTALL = install

CFLAGS = -Wall -I../include -I. #-march=native -mtune=native -ggdb for debug, -O3 for release
LFLAGS = #-march=native -mtune=native

GMPLIB = -L/usr/lib -lgmp 

BIN = bin
OBJ = obj
SRC = src

.PHONY: ensure_dirs clean

all: ensure_dirs $(BIN)/ranked_tree

# Generate object files
$(OBJ)/generate_sarray.o: $(SRC)/generate_sarray.c
	$(CC) $(CFLAGS) -o $@ -c $^

$(OBJ)/hash_table.o: $(SRC)/hash_table.c
	$(CC) $(CFLAGS) -o $@ -c $^

$(OBJ)/utils.o: $(SRC)/utils.c
	$(CC) $(CFLAGS) -o $@ -c $^

$(OBJ)/monitored_memory.o: $(SRC)/monitored_memory.c
	$(CC) $(CFLAGS) -o $@ -c $^

$(OBJ)/newick_tree.o: $(SRC)/newick_tree.c
	$(CC) $(CFLAGS) -o $@ -c $^

$(OBJ)/getopt.o: $(SRC)/getopt.c
	$(CC) $(CFLAGS) -o $@ -c $^

$(OBJ)/lca.o: $(SRC)/lca.c
	$(CC) $(CFLAGS) -o $@ -c $^

# Build binaries
$(BIN)/ranked_tree: $(OBJ)/newick_tree.o $(OBJ)/generate_sarray.o $(OBJ)/utils.o $(OBJ)/getopt.o $(OBJ)/monitored_memory.o $(OBJ)/hash_table.o $(OBJ)/lca.o
	$(LINK) $(LFLAGS) -o $@ $^

ensure_dirs:
	@if [ ! -d $(BIN) ]; then mkdir $(BIN); fi
	@if [ ! -d $(OBJ) ]; then mkdir $(OBJ); fi

clean:
	rm -rf $(OBJ)/*.o $(BIN)/ranked_tree

