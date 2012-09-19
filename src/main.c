#include <stdio.h>
#include <string.h>

int train(int argc, char* argv[]);
int seq2qu(int argc, char* argv[]);

int main(int argc, char *argv[])
{
  if(argc < 2)
  {
    fprintf(stderr, "\nUsage: %s <command> [otions]\n\n", argv[0]);
    fprintf(stderr, "Commands:\n");
    fprintf(stderr, "\ttrain    codon usage and phase training via self-orgnizing approach.\n");
    fprintf(stderr, "\tseq2qu   convert the input sequence file (ACGTTTCGT) into a single quanternary sequence (012333123)\n");
    fprintf(stderr, "\t         the input file can contain multiple sequences, and the output file contains only one sequence, without space (eg. \\n).\n");
    return 1;
  }
  
  if(strcmp(argv[1], "train") == 0) return train(argc - 1, argv + 1);
  else if(strcmp(argv[1], "seq2qu") == 0) return seq2qu(argc - 1, argv + 1);
  else
  {
    fprintf(stderr, "ERROR: unrecognized command '%s'\n", argv[1]);
    return 1;
  }
  
  return 0;
}


