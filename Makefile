# $ make
# $ make install (copy into bin directory)

CC = clang
LINK = clang
INSTALL = install

GMPLIB = -L/usr/local/lib -lgmp -lmpfr

CFLAGS = -Wall -Wunused -I../include -I. #-march=native -mtune=native -ggdb for debug, -O3 for release
LFLAGS = $(GMPLIB) #-march=native -mtune=native

BIN = bin
OBJ = obj
SRC = src

.PHONY: ensure_dirs clean

all: ensure_dirs $(BIN)/ranked_tree

# Generate object files
$(OBJ)/main.o: $(SRC)/main.c
	$(CC) $(CFLAGS) -o $@ -c $^

$(OBJ)/hash_table.o: $(SRC)/hash_table.c
	$(CC) $(CFLAGS) -o $@ -c $^

$(OBJ)/traits.o: $(SRC)/traits.c
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
$(BIN)/ranked_tree: $(OBJ)/newick_tree.o $(OBJ)/main.o $(OBJ)/traits.o $(OBJ)/getopt.o $(OBJ)/monitored_memory.o $(OBJ)/hash_table.o $(OBJ)/lca.o
	$(LINK) $(LFLAGS) -o $@ $^

ensure_dirs:
	@if [ ! -d $(BIN) ]; then mkdir $(BIN); fi
	@if [ ! -d $(OBJ) ]; then mkdir $(OBJ); fi

clean:
	rm -rf $(OBJ)/*.o $(BIN)/ranked_tree

