#ifndef __NEWICKFORM_H__
#define __NEWICKFORM_H__

typedef struct newick_child {
	struct newick_node *node;
	struct newick_child *next;
} newick_child;

typedef struct newick_node {
    int id;
	char *taxon;
	float dist;
	int childNum;
	struct newick_child *child;
	struct newick_node *parent;
} newick_node;

typedef struct newick_bin_node {
    int id;
    float dist;
    char *taxon;

    struct newick_bin_node *parent;
    struct newick_bin_node *left_child;
    struct newick_bin_node *right_child;
} newick_bin_node;

newick_node* parseTree(char *str);
void printTree(newick_node *root);

newick_bin_node* parseBinaryTree(char *str);
void printBinaryTree(newick_bin_node *root);

#endif
