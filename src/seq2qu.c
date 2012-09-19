#include <stdio.h>
#include <zlib.h>
#include "kseq.h"
#include "util.h"


KSEQ_INIT(gzFile, gzread)
int seq2qu(int argc, char* argv[])
{
  gzFile fp;
  FILE*   fq;
  kseq_t* seq;
  int i, l, qu;
  
  if(argc != 3)
  {
    fprintf(stderr, "\nUsage: sega seq2qu <input.[fa, fq][.gz]> <out.qu>\n");
    return 1;
  }
  
  fp = gzopen(argv[1], "r");
  if(fp == NULL)
  {
    fprintf(stderr, "Cannot open %s\n", argv[1]);
    return 1;
  }

  fq = fopen_or_die(argv[2], "w");

  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) 
  {
    for(i = 0; i < l; i++)
    {
      switch(seq->seq.s[i])
      {
      case 'A': qu = 0; break;
      case 'C': qu = 1; break;
      case 'G': qu = 2; break;
      case 'T': qu = 3; break;
      case 'R': qu = (rand() % 2) == 0 ? 0 : 2; break;
      case 'Y': qu = (rand() % 2) == 0 ? 1 : 3; break;
      case 'M': qu = (rand() % 2) == 0 ? 0 : 1; break;
      case 'K': qu = (rand() % 2) == 0 ? 2 : 3; break;
      case 'W': qu = (rand() % 2) == 0 ? 3 : 0; break;
      case 'S': qu = (rand() % 2) == 0 ? 1 : 2; break;
      case 'B': qu = (rand() % 3) == 0 ? 1 : ((rand() % 2) == 0 ? 2 : 3); break;
      case 'D': qu = (rand() % 3) == 0 ? 0 : ((rand() % 2) == 0 ? 2 : 3); break;
      case 'H': qu = (rand() % 3) == 0 ? 1 : ((rand() % 2) == 0 ? 0 : 3); break;
      case 'V': qu = (rand() % 3) == 0 ? 1 : ((rand() % 2) == 0 ? 2 : 1); break;
      case 'N': qu = (rand() % 4) == 0 ? 0 : ((rand() % 3) == 0 ? 1 : (rand() % 2) == 0 ? 2 : 3); break;
      default: qu = 4; suicidef("I cant recognize: %c\n", seq->seq.s[i]); break;
      }
      
      fprintf(fq, "%d", qu);
    }
  }

  kseq_destroy(seq);
  gzclose(fp);
  fclose(fq);
  return 0;
}
