# FK_Simulator


The R files are fairly easy to follow and run. On my laptop it takes 43 seconds to run.

The C code will require a GPU. To compile are run the code, edit the Makefile:
 - edit the path to your CUDA libraries
 - edit the path to your nvcc compiler
 - edit the gpu_architecture flag according to your GPU (gpu_architecture = sm_60 means GPU has compute capability 6.0) see https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#gpu-feature-list
 
For comparison, the C code takes about 0.6 seconds to run.

IMPORTANT:
The C code will not run for ALL combinations of NOBS and NPATHS. This is a limitation of the random number generator in CUDA (and possibly due to my weak coding ability).
 
 I've included MATLAB code for plotting.
