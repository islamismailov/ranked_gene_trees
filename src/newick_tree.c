#include <stdio.h>
#include <string.h>

#include "monitored_memory.h"
#include "newick_tree.h"

newick_node* parseTree(char *str)
{
	newick_node *node;
	newick_child *child;
	char *pcCurrent;
	char *pcStart;
	char *pcColon = NULL;
	char cTemp;
	int iCount;

	pcStart = str;

	if (*pcStart != '(')
	{
		// Leaf node. Separate taxon name from distance. If distance not exist then take care of taxon name only
		pcCurrent = str;
		while (*pcCurrent != '\0')
		{
			if (*pcCurrent == ':')
			{
				pcColon = pcCurrent;
			}
			pcCurrent++;
		}
		node = (newick_node *) monitored_calloc(sizeof(newick_node), 1);
		if (pcColon == NULL)
		{
			// Taxon only
		    node->taxon = (char *) monitored_calloc(sizeof(char), strlen(pcStart) + 1);
			memcpy(node->taxon, pcStart, strlen(pcStart));
		}
		else
		{
			// Taxon
			*pcColon = '\0';
			node->taxon = (char *) monitored_calloc(sizeof(char), strlen(pcStart) + 1);
			memcpy(node->taxon, pcStart, strlen(pcStart));
			*pcColon = ':';
			// Distance
			pcColon++;
			node->dist = (float)atof(pcColon);
		}
		node->childNum = 0;
	}
	else
	{
		// Create node
		node = (newick_node *) monitored_calloc(sizeof(newick_node), 1);
		child = NULL;
		// Search for all child nodes
		// Find all ',' until corresponding ')' is encountered
		iCount = 0;
		pcStart++;
		pcCurrent = pcStart;
		while (iCount >= 0)
		{
			switch (*pcCurrent)
			{
				case '(':
					// Find corresponding ')' by counting
					pcStart = pcCurrent;
					pcCurrent++;
					iCount++;
					while (iCount > 0) // && *pcCurrent != '\0'
					{
						if (*pcCurrent == '(')
						{
							iCount++;
						}
						else if (*pcCurrent == ')')
						{
							iCount--;
						}
						pcCurrent++;
					}
					while (*pcCurrent != ',' && *pcCurrent != ')') // && *pcCurrent != '\0'
					{
						pcCurrent++;
					}
 					cTemp = *pcCurrent;
					*pcCurrent = '\0';
					// Create a child node
					if (child == NULL)
					{
					    node->child = (newick_child *) monitored_calloc(sizeof(newick_child), 1);
						node->childNum = 1;
						child = node->child;
					}
					else
					{
						child->next = (newick_child *) monitored_calloc(sizeof(newick_child), 1);
						node->childNum++;
						child = child->next;
					}
					child->node = parseTree(pcStart);
					*pcCurrent = cTemp;
					if (*pcCurrent != ')')
					{
						pcCurrent++;
					}
				break;

				case ')':
					// End of the tree. Go to next part to retrieve distance
					iCount--;
				break;

				case ',':
					// Impossible separation since according to the algorithm, this symbol will never encountered.
					// Currently don't handle this and don't create any node
                    //throw runtime_error("bad tree");
				break;

				default:
					// leaf node encountered
					pcStart = pcCurrent;
					while (*pcCurrent != ',' && *pcCurrent != ')')
					{
						pcCurrent++;
					}
					cTemp = *pcCurrent;
					*pcCurrent = '\0';
					// Create a child node
					if (child == NULL)
					{
					    node->child = (newick_child *) monitored_calloc(sizeof(newick_child),1);
						node->childNum = 1;
						child = node->child;
					}
					else
					{
					    child->next = (newick_child *) monitored_calloc(sizeof(newick_child),1);
						node->childNum++;
						child = child->next;
					}
					child->node = parseTree(pcStart);
					*pcCurrent = cTemp;
					if (*pcCurrent != ')')
					{
						pcCurrent++;
					}
				break;
			}
		}

		// If start at ':', then the internal node has no name.
		pcCurrent++;
		if (*pcCurrent == ':')
		{
			pcStart = pcCurrent + 1;
			while (*pcCurrent != '\0' && *pcCurrent != ';')
			{
				pcCurrent++;
			}
			cTemp = *pcCurrent;
			*pcCurrent = '\0';
			node->dist = (float)atof(pcStart);
			*pcCurrent = cTemp;
		}
		else if (*pcCurrent != ';' && *pcCurrent != '\0')
		{
			// Find ':' to retrieve distance, if any.
			// At this time *pcCurrent should equal to ')'
			pcStart = pcCurrent;
			while (*pcCurrent != ':')
			{
				pcCurrent++;
			}
			cTemp = *pcCurrent;
			*pcCurrent = '\0';
			node->taxon = (char *) monitored_calloc(sizeof(char), strlen(pcStart) + 1);
			memcpy(node->taxon, pcStart, strlen(pcStart));
			*pcCurrent = cTemp;
			pcCurrent++;
			pcStart = pcCurrent;
			while (*pcCurrent != '\0' && *pcCurrent != ';')
			{
				pcCurrent++;
			}
			cTemp = *pcCurrent;
			*pcCurrent = '\0';
			node->dist = (float)atof(pcStart);
			*pcCurrent = cTemp;
		}
	}

	return node;
}

void printTree(newick_node *root)
{
	newick_child *child;
	if (root->childNum == 0)
	{
		printf("%s:%0.6f", root->taxon, root->dist);
	}
	else
	{
		child = root->child;
		printf("(");
		while (child != NULL)
		{
			printTree(child->node);
			if (child->next != NULL)
			{
				printf(",");
			}
			child = child->next;
		}
		if (root->taxon != NULL)
		{
			printf(")%s:%0.6f", root->taxon, root->dist);
		}
		else
		{
			printf("):%0.6f", root->dist);
		}
	}
}

newick_bin_node* parseBinaryTree(char *str) {
    newick_bin_node *node;

    char *pcCurrent;
    char *pcStart;
    char *pcColon = NULL;
    char cTemp;
    int iCount;

    pcStart = str;

    if (*pcStart != '(')
    {
        // Leaf node. Separate taxon name from distance. If distance not exist then take care of taxon name only
        pcCurrent = str;
        while (*pcCurrent != '\0')
        {
            if (*pcCurrent == ':')
            {
                pcColon = pcCurrent;
            }
            pcCurrent++;
        }
        node = (newick_bin_node *)monitored_calloc(sizeof(newick_bin_node), 1);
        node->parent = node->left_child = node->right_child = NULL;

        if (pcColon == NULL)
        {
            // Taxon only
            node->taxon = (char*)monitored_calloc(sizeof(char), strlen(pcStart) + 1);
            memcpy(node->taxon, pcStart, strlen(pcStart));
        }
        else
        {
            // Taxon
            *pcColon = '\0';
            node->taxon = (char*)monitored_calloc(sizeof(char), strlen(pcStart) + 1);
            memcpy(node->taxon, pcStart, strlen(pcStart));
            *pcColon = ':';
            // Distance
            pcColon++;
            node->dist = (float)atof(pcColon);
        }
    }
    else
    {
        // Create node
        node = (newick_bin_node*)monitored_calloc(sizeof(newick_bin_node), 1);
        node->parent = node->left_child = node->right_child = NULL;

        // Search for all child nodes
        // Find all ',' until corresponding ')' is encountered
        iCount = 0;
        pcStart++;
        pcCurrent = pcStart;
        while (iCount >= 0)
        {
            switch (*pcCurrent)
            {
                case '(':
                    // Find corresponding ')' by counting
                    pcStart = pcCurrent;
                    pcCurrent++;
                    iCount++;
                    while (iCount > 0) // && *pcCurrent != '\0'
                    {
                        if (*pcCurrent == '(')
                        {
                            iCount++;
                        }
                        else if (*pcCurrent == ')')
                        {
                            iCount--;
                        }
                        pcCurrent++;
                    }
                    while (*pcCurrent != ',' && *pcCurrent != ')') // && *pcCurrent != '\0'
                    {
                        pcCurrent++;
                    }
                    cTemp = *pcCurrent;
                    *pcCurrent = '\0';

                    // Parse two first children, ignore everything else
                    if (node->left_child == NULL) {
                        node->left_child = parseBinaryTree(pcStart);
                    } else if (node->right_child == NULL) {
                        node->right_child = parseBinaryTree(pcStart);
                    } else {
                        puts("newick binary tree parser: number of children exceeds 2 for some nodes");
                    }

                    *pcCurrent = cTemp;
                    if (*pcCurrent != ')')
                    {
                        pcCurrent++;
                    }
                break;

                case ')':
                    // End of the tree. Go to next part to retrieve distance
                    iCount--;
                break;

                case ',':
                    // Impossible separation since according to the algorithm, this symbol will never encountered.
                    // Currently don't handle this and don't create any node
                break;

                default:
                    // leaf node encountered
                    pcStart = pcCurrent;
                    while (*pcCurrent != ',' && *pcCurrent != ')')
                    {
                        pcCurrent++;
                    }
                    cTemp = *pcCurrent;
                    *pcCurrent = '\0';

                    // Parse two first children, ignore everything else
                    if (node->left_child == NULL) {
                        node->left_child = parseBinaryTree(pcStart);
                    } else if (node->right_child == NULL) {
                        node->right_child = parseBinaryTree(pcStart);
                    } else {
                        puts("newick binary tree parser: number of children exceeds 2 for some nodes");
                    }

                    *pcCurrent = cTemp;
                    if (*pcCurrent != ')')
                    {
                        pcCurrent++;
                    }
                break;
            }
        }

        // If start at ':', then the internal node has no name.
        pcCurrent++;
        if (*pcCurrent == ':')
        {
            pcStart = pcCurrent + 1;
            while (*pcCurrent != '\0' && *pcCurrent != ';')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            node->dist = (float)atof(pcStart);
            *pcCurrent = cTemp;
        }
        else if (*pcCurrent != ';' && *pcCurrent != '\0')
        {
            // Find ':' to retrieve distance, if any.
            // At this time *pcCurrent should equal to ')'
            pcStart = pcCurrent;
            while (*pcCurrent != ':')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            node->taxon = (char *) monitored_calloc(sizeof(char), strlen(pcStart) + 1);
            memcpy(node->taxon, pcStart, strlen(pcStart));
            *pcCurrent = cTemp;
            pcCurrent++;
            pcStart = pcCurrent;
            while (*pcCurrent != '\0' && *pcCurrent != ';')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            node->dist = (float)atof(pcStart);
            *pcCurrent = cTemp;
        }
    }

    return node;
}

void printBinaryTree(newick_bin_node *root) {
    if (root->left_child == NULL && root->right_child == NULL) {
        printf("%s:%0.6f", root->taxon, root->dist);
    } else {
        printf("(");

        if (root->left_child != NULL) printBinaryTree(root->left_child);
        printf(",");

        if (root->right_child != NULL) printBinaryTree(root->right_child);

        if (root->taxon != NULL) {
            printf(")%s:%0.6f", root->taxon, root->dist);
        } else {
            printf("):%0.6f", root->dist);
        }
    }
}
