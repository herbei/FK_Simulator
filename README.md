# FK_Simulator


The R files are fairly easy to follow and run. On my laptop it takes 43 seconds to run (with NOBS = 200 and NPATHS=200). I am sure the code can be improved a lot -- R is not my forte though.

The C code will require a GPU. To compile are run the code, edit the Makefile:
 - edit the path to your CUDA libraries
 - edit the path to your nvcc compiler
 - edit the gpu_architecture flag according to your GPU (gpu_architecture = sm_60 means GPU has compute capability 6.0) see https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#gpu-feature-list
 
For comparison, the C code takes about 0.6 seconds to run. In the current version, it is difficult / impossible to use the R code for calibration, but the C code will do fine.

IMPORTANT:
The C code will not run for ALL combinations of NOBS and NPATHS. This is a limitation of the random number generator in CUDA (and possibly due to my weak coding ability). I will look into this further, but for now it seems that NOBS = 200 and NPATHS = 200 is the max. One could increase one variable and decrease the other and it the code will likely run OK. Apologies that this is too vague. This is also related to the long-standing question: what is the optimal combination of blocks / threads per block to use - I don't have an answer for that.
 
 I've included MATLAB code for plotting.
