#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <time.h>
#include <curand.h>
#include <curand_kernel.h>



#define NUM_THREADS 1000
#define NUM_BLOCKS 100


#define HH 1e7

typedef struct {
    int width;
    int height;
    double* elements;
} Matrix;


extern "C" void fkpaths(double *domain, Matrix SITES, Matrix OXY, Matrix UV, double *KXY, Matrix FKSOL);

#define CUDA_CALL(x) do { if ((x) != cudaSuccess) { \
printf("Error at %s : %d \n",__FILE__, __LINE__);\
return EXIT_FAILURE;}} while(0)

#define CURAND_CALL(x) do { if ((x) != CURAND_STATUS_SUCCESS) { \
printf("Error at %s : %d\n",__FILE__, __LINE__);\
return EXIT_FAILURE;}} while(0)


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}




//###################################################################
__global__ void setup_kernel(curandState *state)
{


  int bid=blockIdx.x;
  int tid=threadIdx.x;
  int NPATHS = blockDim.x;


  int thread=bid*NPATHS+tid;

  // each thread gets same seed, different seq number, no offset
  curand_init(1234,thread,0,&state[thread]);
}
//##################################################################




__global__ void gpu_fkpaths(curandState *state, double *dev_domain, Matrix dev_OXY, Matrix dev_UV, double *dev_KXY, Matrix dev_SITES, Matrix dev_FKSOL){

	int bid = blockIdx.x;
	int tid = threadIdx.x;
	int NPATHS = blockDim.x;


	int thread = bid*NPATHS + tid;

	// copy state to local memory for efficiency;
  	curandState localState=state[thread];



  	int M,N;
  	M = dev_OXY.height;
  	N = dev_OXY.width;


	double xc,yc; //x,y current
  	double xn,yn; //x,y new
  	double tau=0.0;

  	double delx;
  	delx=(dev_domain[1]-dev_domain[0])/(double)(N-1);
  
  	double dely;
  	dely=(dev_domain[3]-dev_domain[2])/(double)(M-1);


  	xc = dev_SITES.elements[bid*dev_SITES.width + 0];
  	yc = dev_SITES.elements[bid*dev_SITES.width + 1];

  	int i,j;
  	double uc,vc;

  	tau = 0.0;
  	while(  ((xc - dev_domain[0])*(xc - dev_domain[1])<0) && ((yc - dev_domain[2])*(yc - dev_domain[3])<0)  ){

  		//find uvindx
    	j = ceil( (xc - dev_domain[0])/delx );
    	i = ceil( (yc - dev_domain[2])/dely );
    	if (j<0){
      			j=0;
    		}
    	else{
      			if (j>(N-1)) j=N-1;
    		}
    
    	if (i<0){
      			i=0;
    		}
    	else{
      			if (i>(M-1)) i=M-1;
    		}

	    uc=-dev_UV.elements[i*dev_UV.width+j];
    	vc=-dev_UV.elements[i*dev_UV.width+N+j];


    	xn = xc + HH * uc + sqrt(HH)*sqrt(2*dev_KXY[0])*curand_normal(&localState);
    	yn = yc + HH * vc + sqrt(HH)*sqrt(2*dev_KXY[1])*curand_normal(&localState);

    	xc=xn;
    	yc=yn;
    	tau = tau + HH;


  	}


  	double lam = 1e-11;
  	int II,JJ;
  	JJ = ceil( (xc - dev_domain[0])/delx );
  	II = ceil( (yc - dev_domain[2])/dely );

  	if (JJ<0){
    		JJ=0;
  		}
  	else{
    		if (JJ>(N-1)) JJ=N-1;
  		}

  	if (II<0){
    		II=0;
  		}
  	else{
    		if (II>(M-1)) II=M-1;
  		}



  	dev_FKSOL.elements[bid * dev_FKSOL.width + tid] = dev_OXY.elements[II * dev_OXY.width + JJ] * exp(-lam*tau);
	//dev_FKSOL.elements[bid * dev_FKSOL.width + tid] = (double)bid;


  	// copy state back to global memory
  	state[thread]=localState;

}




void fkpaths(double *domain, Matrix SITES, Matrix OXY, Matrix UV, double *KXY, Matrix FKSOL){


	int NOBS;
	int NPATHS;

	NOBS = FKSOL.height;
	NPATHS = FKSOL.width;


 	double *dev_domain;
  	gpuErrchk( cudaMalloc( (void **)&dev_domain, 4*sizeof(double)) );
  	gpuErrchk( cudaMemcpy(dev_domain, domain, 4*sizeof(double), cudaMemcpyHostToDevice) );
  	
	Matrix dev_OXY;
  	dev_OXY.height=OXY.height;
  	dev_OXY.width=OXY.width;
  	gpuErrchk( cudaMalloc( (void **)&dev_OXY.elements, dev_OXY.height*dev_OXY.width*sizeof(double)) );
  	gpuErrchk( cudaMemcpy(dev_OXY.elements, OXY.elements, dev_OXY.height*dev_OXY.width*sizeof(double), cudaMemcpyHostToDevice) );

  	Matrix dev_UV;
	dev_UV.height=UV.height;
  	dev_UV.width=UV.width;
  	gpuErrchk( cudaMalloc( (void **)&dev_UV.elements, dev_UV.height*dev_UV.width*sizeof(double)) );
  	gpuErrchk( cudaMemcpy(dev_UV.elements, UV.elements, dev_UV.height*dev_UV.width*sizeof(double), cudaMemcpyHostToDevice) );

 	double *dev_KXY;
  	gpuErrchk( cudaMalloc( (void **)&dev_KXY, 2*sizeof(double)) );
  	gpuErrchk( cudaMemcpy(dev_KXY, KXY, 2*sizeof(double), cudaMemcpyHostToDevice) );


  	Matrix dev_SITES;
  	dev_SITES.height = SITES.height;
  	dev_SITES.width = SITES.width;
  	gpuErrchk( cudaMalloc( (void **)&dev_SITES.elements, dev_SITES.height*dev_SITES.width*sizeof(double) ) );
  	gpuErrchk( cudaMemcpy(dev_SITES.elements, SITES.elements, dev_SITES.height * dev_SITES.width * sizeof(double), cudaMemcpyHostToDevice) );

	Matrix dev_FKSOL;
  	dev_FKSOL.height=NOBS; 
  	dev_FKSOL.width=NPATHS;
  	gpuErrchk( cudaMalloc( (void **)&dev_FKSOL.elements, dev_FKSOL.height*dev_FKSOL.width*sizeof(double)) );
	gpuErrchk( cudaPeekAtLastError() );
    
  	



  	curandState *devStates;
  	gpuErrchk ( cudaMalloc( (void **)&devStates, NPATHS*sizeof(curandState)) );
  	setup_kernel<<<NOBS, NPATHS>>>(devStates);
	gpuErrchk( cudaPeekAtLastError() );
	

  	gpu_fkpaths<<<NOBS,NPATHS>>>(devStates, dev_domain, dev_OXY, dev_UV, dev_KXY, dev_SITES, dev_FKSOL);
	gpuErrchk( cudaPeekAtLastError() );

	gpuErrchk( cudaMemcpy(FKSOL.elements, dev_FKSOL.elements, FKSOL.height * FKSOL.width * sizeof(double), cudaMemcpyDeviceToHost) );
	gpuErrchk( cudaPeekAtLastError() );
	//gpuErrchk( cudaDeviceSynchronize() );



	//printf("test = %f\n",FKSOL.elements[3]);

  	//free
  	cudaFree(dev_domain);
	cudaFree(dev_OXY.elements);
	cudaFree(dev_UV.elements);
	cudaFree(dev_KXY);
	cudaFree(dev_SITES.elements);
	cudaFree(dev_FKSOL.elements);



	//printf("%s\n", "done.");
}





