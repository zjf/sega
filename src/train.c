#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <time.h>
#include <stdarg.h>
#include "util.h"

static void cnt(FILE *fp, long n, int dn, int dn_read, int cuNum, float p[cuNum+1][64], float **id_score, float *class_prob, float *phase_prob) {
  int i, j, i0, i1, i2, j1, id, id2cu;
  char s[1000];
  float x, id_score_norm[6*cuNum+1];
  
  x = 0.;
  
  for(id = 0; id < 6 * cuNum + 1; id++)
    x += id_score[n/dn][id];
  
  for(id = 0; id < 6 * cuNum + 1; id++) {
    id_score_norm[id] = id_score[n/dn][id] / x;
    phase_prob[id] += id_score_norm[id];
  }
  
  for(id = 0; id < 6 * cuNum + 1; id++) {
    id2cu = ceil(id/6.);    
    class_prob[id2cu] += id_score_norm[id];
    
    j = id%3;
    fseek(fp, n+j, 0);
    fread(s, 1, dn_read, fp);
    for(i = 0; i < dn_read; i++) {
      i0 = s[i++] - '0';
      i1 = s[i++] - '0';
      i2 = s[i] - '0';
      j1 = ((id - (id != 0) * (id2cu-1) * 6)/4 ? 63-i2*16-i1*4-i0 : i0*16+i1*4+i2);
      //printf("i = %d i0 = %d i1 = %d i2 = %d j1 = %d\n", i, i0, i1, i2, j1);
      p[id2cu][j1] += id_score_norm[id];
    }
  }  
}

static float calProb(FILE *fp, long n, int dn, int dn_read, int cuNum, float q[cuNum+1][64], int id) {
  long i, j;
  int i0, i1, i2, j1, id2cu;
  char s[1000];
  float x;

  fseek(fp, n, 0);
  fread(s, sizeof(char), dn, fp);
  
  x = 1.;

  id2cu = ceil(id/6.);

  j = id%3;
  
  for(i=j; i<dn_read+j; i++) {
    i0 = s[i++] -'0';
    i1 = s[i++] - '0';
    i2 = s[i] - '0';
    j1 = ((id - (id != 0) * (id2cu - 1) * 6)/4 ? 63-i2*16-i1*4-i0 : i0*16+i1*4+i2);
    x *= q[id2cu][j1];
  }

  return x;
  
}

static void iter_E(FILE *fp, int n0, int dn, float **id_score, int cuNum, float q[cuNum+1][64], long nw, float *class_prob, float *phase_prob)
{
  long i;
  int j, k, dn_read;
  long n;
  float x, p[cuNum+1][64], p_tc[cuNum+1][64];
  
  //count p, class_prob and phase_prob
  for(i = 0; i < cuNum + 1; i++)
    class_prob[i] = 0.;
  for(i = 0; i< 6 * cuNum + 1; i++)
    phase_prob[i] = 0.;
  for(i = 0; i < cuNum + 1; i++)
    for(j = 0; j < 64; j++)
      p[i][j] = 0;

  n = n0;
  dn_read = dn - 3 + (dn % 3 != 0) * (dn % 3 == 2? 1: -1);
  for(i=0; i<nw; i++) {
    cnt(fp, n, dn, dn_read, cuNum, p, id_score, class_prob, phase_prob);
    n += dn;
    if((i+1)%1000000==0)
      printf("cnt window NO.%ld\n", i);
  }
  //normalize class_prob
  x = 0.;
  for(i=0; i<(cuNum+1); i++)
    x += class_prob[i];
  for(i=0; i<(cuNum+1); i++)
    class_prob[i] = class_prob[i] / x;

  phase_prob[0] = 1;
  for(i=0; i<cuNum; i++) {
    x = 0.;
    for(j=1; j<7; j++)
      x += phase_prob[i*6 + j];
    for(j=1; j<7; j++)
      phase_prob[i*6 + j] = phase_prob[i*6 + j] / x;
  }
  
  //symmetrize non-coding table of p
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      for(k=0;k<4;k++)
        p_tc[0][i*16+j*4+k] = p[0][63-k*16-j*4-i];
  
  
  for(j=0;j<64;j++)
    p[0][j] = (p[0][j] + p_tc[0][j])/2;
  
  
  //convert q to q(order-1)
  for(i=0; i<(cuNum+1); i++) {
    x = 0.;
    for(j=0; j<64; j++)
      x += p[i][j];
    for(j=0; j<64; j++){
      q[i][j] = 64*p[i][j]/x;
    } 
  }
}

static float iter_M(FILE* fp, long n0, int dn, float** id_score, int cuNum, float q[cuNum+1][64], long nw, float* class_prob, float* phase_prob)
{
  int i, j, cuid, id, id2cu, dn_read;
  long n;
  float max, xcu[cuNum], Q = 0.;
  
  dn_read = dn - 3 + (dn % 3 != 0) * (dn % 3 == 2? 1: -1);  
  //calculate posterior probability
  n = n0;
  for(i=0; i<nw; i++) {
    for(j=0; j<(6*cuNum+1); j++) {
      id2cu = ceil(j/6.);
      id_score[i][j] = class_prob[id2cu] * phase_prob[j] * calProb(fp, n, dn, dn_read, cuNum, q, j);
    }
    //Get the max score for calculate Q
    for(j=0; j<(cuNum+1); j++)
      xcu[j] = 0.;
    
    for(j=0; j<(cuNum*6+1); j++) {
      id2cu = ceil(j/6.);
      xcu[id2cu] += id_score[i][j];
    }

    max = xcu[0];
    cuid = 0;
    for(j=1; j<cuNum; j++) {
      if(max<xcu[j]) {
        max = xcu[j];
        cuid = j;
      }
    }

    if(0 == cuid){
      Q += log(max);
    }else{
      max = id_score[i][cuid*6-5];
      id = cuid*6-5;
      for(j=cuid*6-4; j<(1+cuid*6); j++) {
        if(max<id_score[i][j]) {
          max = id_score[i][j];
          id = j;
        }
      }
      Q += log(max);
    }
    
    n += dn;
    if((i+1)%1000000==0)
      printf("infer window NO.%d\n", i);
  }
  return Q;
}

static int err(float e, int cuNum, float q[cuNum+1][64], float q1[cuNum+1][64]) {
  int i, j;
  float x=0.;

  for(i=0; i<(cuNum+1); i++)
    for(j=0; j<64; j++)
      x += (q[i][j] - q1[i][j]) * (q[i][j] - q1[i][j]) / (q[i][j] + q1[i][j] + 0.00000001);
  //x += abs(q[i][j] - q1[i][j]);
  
  for(i=0; i<(cuNum+1); i++)
    for(j=0; j<64; j++)
      q1[i][j] = q[i][j];
  printf("err: %7.4f\n", x);
  return(x>e?1:0);
}

static void fprintf_prob(FILE* fq, float P, float Q, float* class_prob, float* phase_prob, float codon_prob[][64], int cuNum, int k)
{
  int i, j;
  fprintf(fq, "@#### k = %d\n", k);
  fprintf(fq,"P: %7.4f   Q: %7.4f\n", P, Q);
  for(j = 0; j < cuNum + 1; j++)
    fprintf(fq, "%7.4f ", class_prob[j]);
  fprintf(fq, "\n");
    
  fprintf(fq, "%7.4f\n", phase_prob[0]);
  for(i = 0; i < cuNum; i++) {
    for(j = 1; j < 7; j++) {
      fprintf(fq, "%7.4f ", phase_prob[i*6 + j]);
    }
    fprintf(fq, "\n");
  }
      
  for(i = 0; i < cuNum + 1; i++) {
    for(j = 0; j < 64; j++) {
      fprintf(fq, "%7.4f ", codon_prob[i][j]);
    }
    fprintf(fq, "\n");
  }
}

static float calP(int cuNum,float** score, long nw)
{
  long  i;
  int   j;
  float sum;
  float P = 0.;
  
  for(i=0; i<nw; i++) {
    sum = 0.;
    for(j=0; j<(6*cuNum+1); j++) {
      sum += score[i][j];
      // sum += 1;
    }
    P += log(sum);
  }
  return P;
}


static void readProb(const char* fileName, int cuNum, float* class_prob, float* phase_prob, float codon_prob[][64], float codon_prob_old[][64])
{
  // Read frequency table 
  // The first line contains the number of non-coding classes
  // 
  FILE* fp;
  int   i;
  int   j;
  
  if((fp = fopen(fileName, "r")) == NULL) {
    fprintf(stderr, "\nFile %s not found!\n", fileName);
    exit(0);
  }
  for(i = 0; i < cuNum + 1; i++) {
    fscanf(fp, "%f", &class_prob[i]);
  }

  for(i = 0; i < cuNum * 6 + 1; i++) {
    fscanf(fp, "%f", &phase_prob[i]);
  }
  for(i = 0; i < cuNum + 1; i++) {
    for(j = 0; j < 64; j++) {
      fscanf(fp, "%f", &codon_prob[i][j]);
      codon_prob_old[i][j] = codon_prob[i][j];
    }
  }
  fclose(fp);
}


int train(int argc, char *argv[]) {

  if(11 != argc && 13 != argc) {
    fprintf(stderr, "\nUsage: sega train <infile.qu> <outdir> [options]\n");
    fprintf(stderr, "\nOptions:\n");
    fprintf(stderr, "\t-p <string>  optional, use it when you want to begin with specifying a codon usage file instead of generating a random phase for each window.\n");
    fprintf(stderr, "\t-c <int>     categories for coding pattern.\n");
    fprintf(stderr, "\t             Assuming there are 2 different species in the metagenome, you can use -c 2.\n");
    fprintf(stderr, "\t-k <int>     iteration steps, used in absence of -e.\n");
    fprintf(stderr, "\t-e <float>   iteration threshold (the distance of codon usage tables between two iterations), used in absence of -k.\n");
    fprintf(stderr, "\t-n <int>     number of windows for trainning.\n");
    fprintf(stderr, "\t-w <int>     width of window. choose a multiple of 3, eg. 90.\n");
    fprintf(stderr, "\t             Your sequence should be longer than the product of -n and -w options.\n");
    return 0;
  }

  int i, j, k;
  int cuNum;
  int dn;
  long n0;
  long nw;
  int k_or_e = 1;
  int iterNum;
  float iterThres = 0.1;
  int initProbFirst = 0;
  const char* probFileName = "";
  
  for(i = 1; i < argc; i++)
  {
    if(argv[i][0] != '-')
      continue;
    switch(argv[i][1])
    {
    case 'c':
      cuNum = atoi(argv[++i]);
      break;
    case 'n':
      nw = atol(argv[++i]);
      break;
    case 'w':
      dn = atoi(argv[++i]);
      break;
    case 'k':
      iterNum = atoi(argv[++i]);
      k_or_e = 1;
      break;
    case 'e':
      iterThres = atof(argv[++i]);
      k_or_e = 0;
      break;
    case 'p':
      initProbFirst = 1;
      i++;
      probFileName = argv[i];
      break;
    default:
      fprintf(stderr, "Unknown option: %c\n", argv[i][1]);
      break;
    }
  }
      
  float q[cuNum+1][64], q1[cuNum+1][64];
  FILE *fp;  // input sequence
  FILE *fq;  // output
  FILE *fq1; // 
  float class_prob[cuNum+1], phase_prob[cuNum*6+1];

  
  // Alloc memory for each window's scores of (6*cuNum+1) phases
  float **id_score = (float **)malloc(sizeof(float *) * nw);
  for(i=0; i<nw; i++)
    id_score[i] = malloc((6 * cuNum + 1)* sizeof(float));
  
  if(initProbFirst)
  {
    readProb(probFileName, cuNum, class_prob, phase_prob, q, q1);
  }
  else
  {
    srand(time(NULL));
    for(i = 0; i < nw; i++)
      for(j = 0; j < (7*cuNum); j++)
        id_score[i][j] = rand()%10000;
    for(i = 0; i < cuNum + 1; i++)
      for(j = 0; j < 64; j++)
        q1[i][j] = 1.; 
  }

  fp = fopen_or_die(argv[1], "r");

  if(access(argv[2], 0) == -1)
    mkdir(argv[2], 0755);
  char outcut[1024];
  sprintf(outcut, "%s%s", argv[2], "/cut");
  fq = fopen_or_die(outcut, "w");
  char outid[1024];
  sprintf(outid, "%s%s", argv[2], "/id_score");
  fq1 = fopen_or_die(outid, "w");

  k = 0;
  n0 = 0;
  int errAC;
  float P_score = 0., Q_score = 0.;
  do {
    k++;
    if(initProbFirst)
    {
      Q_score = iter_M(fp, n0, dn, id_score, cuNum, q, nw, class_prob, phase_prob);
      iter_E(fp, n0, dn, id_score, cuNum, q, nw, class_prob, phase_prob);
    }
    else
    {
      iter_E(fp, n0, dn, id_score, cuNum, q, nw, class_prob, phase_prob);
      Q_score = iter_M(fp, n0, dn, id_score, cuNum, q, nw, class_prob, phase_prob);
    }
    
    P_score = calP(cuNum, id_score, nw);
    fprintf_prob(fq, P_score, Q_score, class_prob, phase_prob, q, cuNum, k);
    errAC = err(iterThres, cuNum, q, q1);

  }while(k_or_e ? (k < iterNum) : errAC);
  
  fclose(fp);
  
  /********************************************************
   *  Assign phase with maximal probability to each window
   *********************************************************/  
  int cuid, id, id2cu;
  float max, xcu[cuNum];
  for(j=0; j<nw; j++) {
    
    for(k=0; k<cuNum+1; k++)
      xcu[k] = 0.;

    for(k=0; k<(cuNum*6+1); k++){
      id2cu = ceil(k/6.);
      xcu[id2cu] += id_score[j][k];
    }

    max = xcu[0];
    cuid = 0;
    for(i=1; i<(cuNum+1); i++) {
      if(max<xcu[i]) {
        max = xcu[i];
        cuid = i;
      }
    }
  
    if(cuid == 0) {
      id = 0;
      fprintf(fq1, "%d\t%e\n", id, max);
    }else{
      max = id_score[j][cuid*6-5];
      id = cuid*6-5;
      for(i=cuid*6-4; i<(cuid*6+1); i++) {
        if(max<id_score[j][i]) {
          max = id_score[j][i];
          id = i;
        }
      }
      fprintf(fq1, "%d\t%e\n", id, max);
    }
    
  }
  fclose(fq1);
  

  for(i = 0; i < nw; i++)
    free(id_score[i]);
  free(id_score);

  return 0;
}
