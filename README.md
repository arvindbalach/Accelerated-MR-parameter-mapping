
# Structured matrix completion algorithm for accelerated Magnetic Resonance (MR) parameter mapping
# PUBLICATION
A. Balachandrasekaran, V. Magnotta and M. Jacob. "Recovery of damped exponentials using structured low rank matrix completion." *IEEE Transactions on Medical Imaging*, 36(10):2087-2098, 2017. 
# Codes

The main files inside the proposed_fast and direct_implementation folders solve the following optimization problem

![](https://latex.codecogs.com/gif.latex?%24%5Cmin_%7B%5Cwidehat%20%5Crho%7D%20%5C%7CT%28%5Cwidehat%20%5Crho%29%5C%7C_p%20&plus;%20%28%5Cmu/2%29%5C%7CA%28%5Cwidehat%20%5Crho%29-b%5C%7C%5E2_2%24)

with and without introducing approximations respectively. ![](https://latex.codecogs.com/gif.latex?%24T%28%5Cwidehat%20%5Crho%29%24) is the Toeplitz matrix formed from the Fourier samples ![](https://latex.codecogs.com/gif.latex?%24%5Cwidehat%20%5Crho%24), ![](https://latex.codecogs.com/gif.latex?%24%5C%7CX%5C%7C_p%24) is the Schatten p-norm of ![](https://latex.codecogs.com/gif.latex?%24X%24), ![](https://latex.codecogs.com/gif.latex?%24A%24) is the Fourier undersampling operator and ![](https://latex.codecogs.com/gif.latex?%24b%24) represents undersampled Fourier measurements. Please refer to the above mentioned publication to get more information on the notations.

## Proposed fast codes
### main_coilcombined
This file recovers the coil combined Fourier data from 30% uniform random Fourier measurements. Run this code to replicate the results in Figure 4 of the paper.

### main_multichannel
This file recovers the multi-channel Fourier data from 12 fold undersampled Fourier measurements. Run this code to replicate the results in Figure 5 of the paper.

### initial_epsilon
Code to initialize the parameter ![](https://latex.codecogs.com/gif.latex?%24%5Cepsilon%24) for the IRLS optimization algorithm.

### estimate_mu
Code to estimate the polynomial ![](https://latex.codecogs.com/gif.latex?%24%5Cmu%24), which is to be used to form the matrix D. Note that this matrix appears in equation (26) of the paper.

### pre_computeChC
Code to compute G that appears in equation (26) of the paper.

### regul_hybrid_precompute
code for computing the first term of equation (26) in the paper.

### xsub_hybrid_sc and xsub_hybrid_mc
Codes for solving equation (26) of the paper. They correspond to single and multi-channel cases respectively.

### compute_cost
Code for evaluating the cost function (mentioned above).
