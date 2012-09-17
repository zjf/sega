/*
 *=========================================================================
 *
 *  Filename: scanf_meta_EM.c
 *  Description: 
 *
 *
 *  Version: 1.0
 *  Created: 2011-09-06
 *  Compiler: gcc scanf_meta_EM.c -o scanf_meta_EM.o -lm
 *
 *  Author: zhujf, zhujianfeng@genomics.cn
 *  Orgnization: ST_BAR
 *
 *
 *=========================================================================
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>


void cnt(FILE *fp, long n, int dn, int cuNum, float p[cuNum+1][64], float **id_score, float *class_prior_prob, float *phase_prior_prob) {
  int i, j, i0, i1, i2, j0, j1, id, id2cu;
  char s[1000];
  float x, id_score_norm[6*cuNum+1];
  
  x = 0.;
  
  for(id = 0; id < (6*cuNum+1); id++)
    x += id_score[n/dn][id];
  
  for(id = 0; id < (6*cuNum+1); id++) {
    id_score_norm[id] = id_score[n/dn][id] / x;
    phase_prior_prob[id] += id_score_norm[id];
  }
  

  for(id=0; id < (6*cuNum+1); id++) {

    id2cu = ceil(id/6.);    
    class_prior_prob[id2cu] += id_score_norm[id];
    
    j = id%3;
    fseek(fp, n+j, 0);
    fread(s, 1, dn-3, fp);
    for(i=0; i<dn-3; i++) {
      i0 = s[i++] - '0';
      i1 = s[i++] - '0';
      i2 = s[i] - '0';
      j1 = ((id - (id != 0) * (id2cu-1) * 6)/4 ? 63-i2*16-i1*4-i0 : i0*16+i1*4+i2);
      p[id2cu][j1] += id_score_norm[id];
    //printf("%ld ", p[j0][j1]);
    }
  }
  
  //printf("\n");
  
}

float calProb(FILE *fp, long n, int dn, int cuNum, float q[cuNum+1][64], int id) {
  int i, j, i0, i1, i2, j0, j1, id2cu;
  char s[1000];
  float x;

  fseek(fp, n, 0);
  fread(s, sizeof(char), dn, fp);
  
  x = 1.;

  id2cu = ceil(id/6.);

  j = id%3;
  
  for(i=j; i<dn+j-3; i++) {
    i0 = s[i++] -'0';
    i1 = s[i++] - '0';
    i2 = s[i] - '0';
    j1 = ((id - (id != 0) * (id2cu - 1) * 6)/4 ? 63-i2*16-i1*4-i0 : i0*16+i1*4+i2);
    x *= q[id2cu][j1];
  }

  return x;
  
}


float iter(FILE *fp, int n0, int dn, float **id_score, int cuNum, float q[cuNum+1][64], int nw, float *class_prior_prob, float *phase_prior_prob) {
//n0: init site, nw: # of win, dn: win width

  int i, j, k, l, cuid, id, id2cu;
  long n;
  float x, p[cuNum+1][64], p_tc[cuNum+1][64], max, xcu[cuNum+1], Q=0.;

  //not necessary
  for(i=0; i<(cuNum+1); i++)
    for(j=0; j<64; j++)
      p[i][j] = rand()%10+1.;
 
  //debug: keep class_prior_prob unchanged
  //  for(i=0; i<cuNum; i++)
  //    class_prior_prob[i] = 1./3;
  
 
  //calculate posterior probability
  n = n0;
  for(i=0; i<nw; i++) {
    for(j=0; j<(6*cuNum+1); j++) {
      id2cu = ceil(j/6.);
      id_score[i][j] = class_prior_prob[id2cu] * phase_prior_prob[j] * calProb(fp, n, dn, cuNum, q, j);
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

 //count p, class_prior_prob and phase_prior_prob
  for(i=0; i<(cuNum+1); i++)
    class_prior_prob[i] == 0.;
  for(i=0; i<(6*cuNum+1); i++)
    phase_prior_prob[i] == 0.;
  
  n = n0;
  for(i=0; i<nw; i++) {
    cnt(fp, n, dn, cuNum, p, id_score, class_prior_prob, phase_prior_prob);
    n += dn;
    if((i+1)%1000000==0)
      printf("cnt window NO.%d\n", i);
  }

  //normalize class_prior_prob
  x = 0.;
  for(i=0; i<(cuNum+1); i++)
    x += class_prior_prob[i];
  for(i=0; i<(cuNum+1); i++)
    class_prior_prob[i] = class_prior_prob[i] / x;

  phase_prior_prob[0] = 1;
  for(i=0; i<cuNum; i++) {
    x = 0.;
    for(j=1; j<7; j++)
      x += phase_prior_prob[i*6 + j];
    for(j=1; j<7; j++)
      phase_prior_prob[i*6 + j] = phase_prior_prob[i*6 + j] / x;
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

  return Q;
  
  
}



int err(float e, int cuNum, float q[cuNum+1][64], float q1[cuNum+1][64]) {
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


int main(int argc, char *argv[]) {

  if(10 != argc) {
    fprintf(stderr, "\nUsage: ./scanf_meta_iniCut_EM.o [infile.qu] [infile.Pr_cut] [outdir] [-c codon_usage_number] [-k/e (# of iter)/(convergence criteria)] [-n # of window]\n");
    return 0;
  }

  int i, i1, j, k, m, dn, n0, cuNum, nw;
  cuNum = atoi(argv[5]);
  nw = atoi(argv[9]);
  float q[cuNum+1][64], q1[cuNum+1][64];
  FILE *fp, *fp1, *fq, *fq1;
  float class_prior_prob[cuNum+1], phase_prior_prob[cuNum*6+1];

  phase_prior_prob[0] = 1;
  for(i=1; i<(6*cuNum+1); i++)
    phase_prior_prob[i] = 1./6;
  

  float **id_score = (float **)malloc(sizeof(float *) * nw);
  for(i=0; i<nw; i++)
    id_score[i] = malloc((6 * cuNum + 1)* sizeof(float));
  
  
  dn = 90;
  n0 = 0;

  
  if((fp=fopen(argv[1], "r")) == NULL) {
    fprintf(stderr, "\nFile not found!\n");
    exit(0);
  }
  /////////////////////////////////////////////////////////
  if((fp1=fopen(argv[2], "r")) == NULL) {
    fprintf(stderr, "\nFile not found!\n");
    exit(0);
  }
  for(i=0; i < cuNum+1; i++) {
    fscanf(fp1, "%f", &class_prior_prob[i]);
    printf("%f ", class_prior_prob[i]);
  }
  printf("\n");

  for(i=0; i<(cuNum*6+1); i++) {
    fscanf(fp1, "%f", &phase_prior_prob[i]);
    printf("%f ", phase_prior_prob[i]);
  }
  printf("\n");
  for(i=0; i < cuNum+1; i++) {
    for(j=0; j < 64; j++) {
      fscanf(fp1, "%f", &q[i][j]);
      q1[i][j] = q[i][j];
      printf("%f ", q[i][j]);
    }
    printf("\n");
  }
  fclose(fp1);
  //////////////////////////////////////////////////////////
  if(access(argv[3], 0)==-1)
    mkdir(argv[3], 0755);

  char outcut[200];
  sprintf(outcut, "%s%s", argv[3], "/cut");
  fq = fopen(outcut, "w");

  k = 0;
  int k_or_e = (*++argv[6] == 'k') ? 1 : 0;
  int iterNum = atoi(argv[7]);
  int errAC;
  float thres = atof(argv[7]);
  float P_score = 0., Q_score = 0., sum;
  
  do {
    k++;
    Q_score = iter(fp, n0, dn, id_score, cuNum, q, nw, class_prior_prob, phase_prior_prob);

    P_score = 0.;
    for(i=0; i<nw; i++) {
      sum = 0.;
      for(j=0; j<(6*cuNum+1); j++) {
        sum += id_score[i][j];
      }
      P_score += log(sum);
    }

    fprintf(fq,"P: %7.4f   Q: %7.4f\n", P_score, Q_score);
    for(j=0; j<cuNum+1; j++)
      fprintf(fq, "%7.4f ", class_prior_prob[j]);
    fprintf(fq, "\n");
    
    fprintf(fq, "%7.4f\n", phase_prior_prob[0]);
    for(i=0; i<cuNum; i++) {
      for(j=1; j<7; j++) {
        fprintf(fq, "%7.4f ", phase_prior_prob[i*6 + j]);
      }
      fprintf(fq, "\n");
    }
    
      
    for(i=0; i<(cuNum+1); i++) {
      fprintf(fq, "\n");
      for(j=0; j<64; j++) {
        fprintf(fq, "%7.4f ", q[i][j]);
      }
    }
    fprintf(fq, "\n@#### k = %d\n !\n", k);

    errAC = err(thres, cuNum, q, q1);

  }while(k_or_e ? (k < iterNum) : errAC);
  
  fclose(fp);
  
  
  char outid[100];
  sprintf(outid, "%s%s", argv[3], "/id_score");
  fq1 = fopen(outid, "w");
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
  

  for(i=0; i<nw; i++)
    free((void *)id_score[i]);
  free((void *)id_score);
  
  
  return 0;
}
