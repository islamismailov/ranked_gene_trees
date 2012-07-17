/*
 * lca.h
 *
 *  Created on: May 30, 2012
 *      Author: islamismailov
 */

#ifndef LCA_H_
#define LCA_H_

#include "generate_sarray.h"

void lca_init(int n, node2int_array *);
int lca(int a, int b);
void lca_end();

#endif /* LCA_H_ */
