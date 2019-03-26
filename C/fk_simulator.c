#include <stddef.h>  // NULL, size_t
#include <math.h> // expf
#include <stdio.h> // printf
#include <time.h> // time
#include <sys/time.h> // gettimeofday
#include <assert.h>
#include <stdlib.h>



typedef struct {
    int width;
    int height;
    double* elements;
} Matrix;


extern void fkpaths(double *domain, Matrix SITES, Matrix OXY, Matrix UV, double *KXY,  Matrix FKSOL);

#define PI 3.14159265








int main(void){
  

  int NOBS = 200;
  int NPATHS = 200;


  int i,j,k;
  FILE *filein, *outfile;
  double value;

  double xy[4]={-34.0, 11.0, -32.0, -3.0};

  double lat0=-17.5;
  double *domain=(double *)malloc(4*sizeof(double));
  domain[0]=xy[0]*111e3*cos(lat0*PI/180.0);
  domain[1]=xy[1]*111e3*cos(lat0*PI/180.0);
  domain[2]=xy[2]*111e3;
  domain[3]=xy[3]*111e3;


  int N = 37;
  int M = 19;

  double delx;
  delx=(domain[1]-domain[0])/(double)(N-1);

  double dely;
  dely=(domain[3]-domain[2])/(double)(M-1);

  //printf("delx=%f dely=%f\n",delx,dely);



  // read velocity vectors u,v
  Matrix UV;
  UV.height=M;
  UV.width=2*N;
  UV.elements=(double *)malloc(UV.height*UV.width*sizeof(double));
  
  filein=fopen("gridu.txt", "r");
  for(i=0;i<M;i++){
    for(j=0;j<N;j++){
      fscanf(filein, "%lf",&value);
      UV.elements[i*UV.width+j]=value;
    }
  }
  fclose(filein);
  
  filein=fopen("gridv.txt", "r");
  for(i=0;i<M;i++){
    for(j=0;j<N;j++){
      fscanf(filein, "%lf",&value);
      UV.elements[i*UV.width+N+j]=value;
    }
  }
  fclose(filein);
  //------------


  // diffusion coefficients
  double *KXY=(double *)malloc(2*sizeof(double));
  KXY[0]=700.0;
  KXY[1]=200.0;
  //------------



  // grid Xp, Yp
  Matrix Xp;
  Xp.width=N;
  Xp.height=M;
  Xp.elements=(double *)malloc(Xp.height*Xp.width*sizeof(double));

  Matrix Yp;
  Yp.width=N;
  Yp.height=M;
  Yp.elements=(double *)malloc(Yp.height*Yp.width*sizeof(double));

  
  for(i=0;i<M;i++){
    for(j=0;j<N;j++){
      Xp.elements[i*Xp.width+j] = domain[0]+delx*((double)j);
      Yp.elements[i*Yp.width+j] = domain[2]+dely*((double)i);
    }
  }


  // oxygen field (only boundary values are needed, read a whole matrix anyway)
  Matrix OXY;
  OXY.height=M;
  OXY.width=N;
  OXY.elements=(double *)malloc(OXY.height*OXY.width*sizeof(double));

  filein=fopen("oxygengrid.txt", "r");
  for(i=0;i<OXY.height;i++){
    for(j=0;j<OXY.width;j++){
      fscanf(filein, "%lf",&value);
      OXY.elements[i*OXY.width+j]=value;
    }
  }
  fclose(filein);
  //-------------------






  Matrix SITES;
  SITES.height=NOBS;
  SITES.width=2;
  SITES.elements=(double *)malloc(SITES.height*SITES.width*sizeof(double));

  filein=fopen("sites_fk.txt", "r");
  for(k=0;k<NOBS;k++){
    fscanf(filein, "%lf",&value);
    SITES.elements[k*SITES.width + 0] = value;

    fscanf(filein, "%lf",&value);
    SITES.elements[k*SITES.width + 1] = value;

  }


  
  Matrix FKSOL;
  FKSOL.height = NOBS;
  FKSOL.width = NPATHS;
  FKSOL.elements=(double *)malloc(FKSOL.height*FKSOL.width*sizeof(double));

  clock_t start, end;
  double cpu_time_used;
  start = clock();
  fkpaths(domain, SITES, OXY, UV, KXY, FKSOL);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("elapsed time = %f\n",cpu_time_used);




  outfile=fopen("oxygen_fk_gpu.txt", "w");
  for(k=0;k<NOBS;k++){
    for(i=0;i<NPATHS;i++){
      fprintf(outfile, "%f ", FKSOL.elements[k*FKSOL.width + i]);
    }
    fprintf(outfile,"\n");
  }
  fclose(outfile);


  return 1;
}





















