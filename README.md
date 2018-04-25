
# Structured matrix completion algorithm for accelerated Magnetic Resonance (MR) parameter mapping
# PUBLICATION
A. Balachandrasekaran, V. Magnotta and M. Jacob. "Recovery of damped exponentials using structured low rank matrix completion." *IEEE Transactions on Medical Imaging*, 36(10):2087-2098, 2017. 
# Codes

The main files inside the proposed_fast and direct_implementation folders solve the following optimization problem

![](https://latex.codecogs.com/gif.latex?%24%5Cmin_%7B%5Cwidehat%20%5Crho%7D%20%5C%7CT%28%5Cwidehat%20%5Crho%29%5C%7C_p%20&plus;%20%28%5Cmu/2%29%5C%7CA%28%5Cwidehat%20%5Crho%29-b%5C%7C%5E2_2%24)

with and without introducing approximations respectively. ![](https://latex.codecogs.com/gif.latex?%24%5Cmu%24) is the regularization parameter, ![](https://latex.codecogs.com/gif.latex?%24T%28%5Cwidehat%20%5Crho%29%24) is the Toeplitz matrix formed from the Fourier samples ![](https://latex.codecogs.com/gif.latex?%24%5Cwidehat%20%5Crho%24), ![](https://latex.codecogs.com/gif.latex?%24%5C%7CX%5C%7C_p%24) is the Schatten p-norm of ![](https://latex.codecogs.com/gif.latex?%24X%24), ![](https://latex.codecogs.com/gif.latex?%24A%24) is the Fourier undersampling operator and ![](https://latex.codecogs.com/gif.latex?%24b%24) represents undersampled Fourier measurements. Please refer to the above mentioned publication to get more information on the notations.

## Proposed fast codes
### main_coilcombined.m
This file recovers the coil combined Fourier data from 30% uniform random Fourier measurements. Run this code to replicate the results in Figure 4 of the paper.

### main_multichannel.m
This file recovers the multi-channel Fourier data from 12 fold undersampled Fourier measurements. Run this code to replicate the results in Figure 5 of the paper.

### initial_epsilon.m
Code to initialize the parameter ![](https://latex.codecogs.com/gif.latex?%24%5Cepsilon%24) for the IRLS optimization algorithm.

### estimate_mu.m
Code to estimate the polynomial ![](https://latex.codecogs.com/gif.latex?%24%5Cmu%24), which is used to form the matrix D. Note that this matrix appears in equation (26) of the paper.

### pre_computeChC.m
Code to compute G that appears in equation (26) of the paper.

### regul_hybrid_precompute.m
code for computing the first term of equation (26) in the paper.

### xsub_hybrid_sc.m and xsub_hybrid_mc.m
Codes for solving equation (26) of the paper. They correspond to single and multi-channel cases respectively.

### compute_cost.m
Code for evaluating the cost function (mentioned above).

### t2_MapsCalc.m
Estimates the T2 map using the time series of images. 

## Direct implementation codes
### main_coilcombined.m
This file recovers the coil combined Fourier data from 30% uniform random Fourier measurements, without introducing approximations. The results obtained are shown in Figure 3 of the paper.

### initial_epsilon.m
Code to initialize the parameter ![](https://latex.codecogs.com/gif.latex?%24%5Cepsilon%24) for the IRLS optimization algorithm.

### estimateQ.m
Code for estimating the matrix Q, which is then used to solve equation (17) in the paper. For details, refer to the Weight update section in the paper.

### xsub_sc_direct.m
Code for solving equation (17) of the paper. It corresponds to the single channel case.

### regul_direct.m
Computes the gradient of the first term in equation (19) of the paper.

### compute_cost.m
Code for evaluating the cost function (mentioned above).

## Associated codes
### operators_sc(fwd,bwd) and operators_mc(fwd,bwd)
These folders contain codes for defining the forward and backward operators for single and multi-channel cases respectively.

### getkspaceindices
The codes inside the folder return the indices corresponding to the Fourier data and the filter coefficients. The indices returned are consistent with matlab convention.

### direct-liftandunlift-codes
#### im2colstep.m
Code used to form the Toeplitz matrix. Rearranges blocks of data into columns of the matrix.

#### col2imstep.m
Code used to extract the data from the Toeplitz matrix. Rearranges columns of the matrix into blocks.



