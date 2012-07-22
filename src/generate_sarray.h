//
//  generate_sarray.h
//  ranked_gene_tree
//
//  Created by Islam Ismailov on 2012/17/7.
//  Copyright (c) 2012 Islam Ismailov. All rights reserved.
//

#ifndef ranked_gene_tree_generate_sarray_h
#define ranked_gene_tree_generate_sarray_h

#include "utils.h"
#include "newick_tree.h"
#include <stdlib.h>

typedef struct matidx {
    int i, j;
} matidx;

typedef struct node2float {
    newick_node *node;
    float val;
} node2float;

typedef struct node2int {
    newick_node *node;
    int val;
} node2int;

typedef char * char_ptr;
typedef int * int_ptr;
typedef int ** int_ptr_ptr;
typedef int *** int_ptr_ptr_ptr;
typedef newick_node * newick_node_ptr;

DEF_ARRAY_DECL(int);
DEF_ARRAY_DECL(char);
DEF_ARRAY_DECL(float);

DEF_ARRAY_DECL(newick_node_ptr);
DEF_ARRAY_DECL(newick_node_ptr_array);

DEF_ARRAY_DECL(node2float);
DEF_ARRAY_DECL(node2float_array);
DEF_ARRAY_DECL(node2int);
DEF_ARRAY_DECL(char_ptr);
DEF_ARRAY_DECL(int_array);
DEF_ARRAY_DECL(int_array_array);

#endif
