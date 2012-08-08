/*
 * lca.c
 *
 *  Created on: May 30, 2012
 *      Author: islamismailov
 */

#include "lca.h"
#include "newick_tree.h"
#include <stdio.h>
#include <stdlib.h>

int timer, n, l;
int *tin, *tout;
int **up;

void lca_preprocess(node2int_array *coalescence_array, int v, int p) {
    int i, to;
    node2int *cp;
    newick_child *np;
    
    tin[v] = ++timer;
    up[v][0] = p;
    for (i = 1; i <= l; ++i)
        up[v][i] = up[up[v][i - 1]][i - 1];
    
    // iterate through children of coalescence u[v]:
    
    for (np = coalescence_array->array[v].node->child; np != NULL; np = np->next) {
#ifndef NDEBUG
        printf("\titerating thru children of u[%d]@%p:\n", v, coalescence_array->array[v].node);
#endif
        // find child's index (to)
        for (cp = coalescence_array->array, to = 0; cp != coalescence_array->last; ++to, ++cp) {
            //if (cp->node == coalescence_array->array[v].node)
            if (cp->node == np->node)
                break;
        }
#ifndef NDEBUG
        printf("\t\tu[%d](u[%d])@%p is a child? ", to, cp->val, np->node);
#endif
        if (to != p && to != coalescence_array->last - coalescence_array->array) {
#ifndef NDEBUG
            puts("y");
#endif
            lca_preprocess(coalescence_array, cp->val, v);
        } else {
#ifndef NDEBUG
            puts("n");
#endif
        }
    }
    
    tout[v] = ++timer;
}

int upper(int a, int b) {
    return tin[a] <= tin[b] && tout[a] >= tout[b];
}

int lca(int a, int b) {
    int i;
    if (upper(a, b)) return a;
    if (upper(b, a)) return b;
    
    for (i = l; i >= 0; --i) {
        if (!upper(up[a][i], b)) {
            a = up[a][i];
        }
    }
    
    return up[a][0];
}

void lca_init(int n, node2int_array *coalescence_array) {
    int i;
    
    tin = (int *) malloc(n * sizeof(int));
    tout = (int *) malloc(n * sizeof(int));
    up = (int **) malloc(n * sizeof(int *));
    
    l = 1;
    while ((1 << l) <= n)  ++l;
    for (i = 0; i < n; ++i)  up[i] = (int *) malloc((l + 1) * sizeof(int));
    
    lca_preprocess(coalescence_array, 0, 0);
}

void lca_end() {
    int i;
    free(tin);
    free(tout);
    for (i = 0; i < n; ++i) free(up[i]);
    free(up);
}
