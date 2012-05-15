#include "seqUtil.h"
#include "Newickform.h"

#define STR_OUT	"out"

int main(int argc, char *argv[])
{
	int iLen, iMaxLen;
	char *pcTreeStr;
	char *pcInputFile;
	char *pcOutputFile;
	newick_node *root;
	char acStrArray[256];

	FILE *f;
	int i;

	// Initialize memory management procedure
	seqMemInit();

	if (argc < 3)
	{
		printf("Usage: SeqAllocator.exe -i (newick tree file) [-o outputfile]\nNote: If outputfile is not specified, the filename for output will be your input_file.out\n");
		exit(-1);
	}
	pcInputFile = NULL;
	pcOutputFile = NULL;
	sDEBUG = 0;
	for (i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-i") == 0)
		{
			i++;
			pcInputFile = argv[i];
			if (pcOutputFile != NULL)
			{
				seqFree(pcOutputFile);
			}
			pcOutputFile = seqMalloc(strlen(pcInputFile) + strlen(STR_OUT) + 2);
			memset(pcOutputFile, '\0', strlen(pcInputFile) + strlen(STR_OUT) + 2);
			sprintf(pcOutputFile, "%s.%s", argv[i], STR_OUT);
		}
		else if (strcmp(argv[i], "-j") == 0)
		{
			i++;
			if (pcOutputFile != NULL)
			{
				seqFree(pcOutputFile);
			}
			pcOutputFile = argv[i];
		}
		else if (strcmp(argv[i], "-DEBUG") == 0)
		{
			sDEBUG = 1;
		}
		else
		{
			printf("Unrecognized identifier: %s\nProgram stop.\n, argv[i]");
			seqFreeAll();
			exit(-1);
		}
	}
	if (pcInputFile == NULL)
	{
		printf("Please specify input tree file using -i parameter.\n");
		seqFreeAll();
		exit(-1);
	}

	// Open tree file
	f = fopen(pcInputFile, "r+");
	if (f == NULL)
	{
		printf("Cannot open input file. Please check file name again.\n");
		seqFreeAll();
		exit(-1);
	}
	// Read in tree string
	pcTreeStr = NULL;
	iLen = 0;
	iMaxLen = 0;
	while (1)
	{
		memset(acStrArray, '\0', 256);
		fgets(acStrArray, 255, f);
		if (acStrArray[0] == '\0' && feof(f))
		{
			break;
		}
		inputString(acStrArray, &pcTreeStr, &iLen, &iMaxLen);
	}
	fclose(f);

	// Parse tree string
	root = parseTree(pcTreeStr);
	printTree(root);
	printf("\n");

	// Free occupied memory
	seqFree(pcOutputFile);
	seqFree(pcTreeStr);

	// End memory management procedure and free all allocated space
	seqFreeAll();
}

