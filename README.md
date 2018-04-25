
# Structured matrix completion algorithm for accelerated Magnetic Resonance (MR) parameter mapping
# PUBLICATION
A. Balachandrasekaran, V. Magnotta and M. Jacob. "Recovery of damped exponentials using structured low rank matrix completion." *IEEE Transactions on Medical Imaging*, 36(10):2087-2098, 2017. 


# Codes

The main files inside the proposed_fast and direct_implementation folders solve the following optimization problem

![](https://latex.codecogs.com/gif.latex?%24%5Cmin_%7B%5Cwidehat%20%5Crho%7D%20%5C%7CT%28%5Cwidehat%20%5Crho%29%5C%7C_p%20&plus;%20%28%5Cmu/2%29%5C%7CA%28%5Cwidehat%20%5Crho%29-b%5C%7C%5E2_2%24)

with and without introducing approximations respectively. ![](https://latex.codecogs.com/gif.latex?%24%5Cmu%24) is the regularization parameter, ![](https://latex.codecogs.com/gif.latex?%24T%28%5Cwidehat%20%5Crho%29%24) is the Toeplitz matrix formed from the Fourier samples ![](https://latex.codecogs.com/gif.latex?%24%5Cwidehat%20%5Crho%24), ![](https://latex.codecogs.com/gif.latex?%24%5C%7CX%5C%7C_p%24) is the Schatten p-norm of ![](https://latex.codecogs.com/gif.latex?%24X%24), ![](https://latex.codecogs.com/gif.latex?%24A%24) is the Fourier undersampling operator and ![](https://latex.codecogs.com/gif.latex?%24b%24) represents undersampled Fourier measurements. Please refer to the above mentioned publication to get more information on the notations.

**To replicate the results corresponding to the proposed fast method in the paper, run the codes main_coilcombined.m and main_multichannel.m present inside the proposed_fast folder. To see the results corresponding to the direct implementation run the main_coilcombined.m present inside the direct_implementation folder. Both the proposed method and the direct implementation requires additional codes present in the folders  Associated codes and  direct-liftandunlift-codes.**

## Proposed fast codes

The codes basically solve equation (17) and (18) in a very efficient manner. Fourier domain approximations are introduced which eliminate the need to store the Toeplitz matrix.

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
The codes solve equations (17) and (18) directly without introducing any approximations. 

### main_coilcombined.m
This file recovers the coil combined Fourier data from 30% uniform random Fourier measurements, without introducing approximations. The results obtained are shown in Figure 3 of the paper.

### initial_epsilon.m
Code to initialize the parameter ![](https://latex.codecogs.com/gif.latex?%24%5Cepsilon%24) for the IRLS optimization algorithm.

### estimateQ.m
Code for estimating the weight matrix which is needed to solve equation (17) in the paper. For details, refer to the Weight update section in the paper.

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
 Code to rearrange the columns of the Toeplitz matrix into blocks of data. This is used to implement the adjoint operator that appears when the graient of equation (17) is computed.

# Data
A fully sampled axial 2-D dataset was acquired on a Siemens 3T Trio scanner with 12 coils using a turbo spin echo sequence. The scan parameters were: TR = 2500 ms, slice thickness = 5 mm, Matrix size = 128x128 and FOV = 22x22 ![](https://latex.codecogs.com/gif.latex?%24%5Cmbox%7Bcm%7D%5E%7B2%7D%24). By varying the echo times (TE) from 10 to 120 ms, we acquired T2 weighted images at twelve equispaced TE.

For the multi-channel case, the data was retrospectively undersampled using a combination of uniform Cartesian and a pseudo-random variable density sampling patterns. Specifically, we uniformly undersampled the x and y directions by a factor of 2 and refer this sampling mask as a 2x2 uniform Cartesian mask. To increase the incoherence between the frames, we also shifted every frame of the mask by zero or one unit (done randomly) along the x and y directions. We achieved an acceleration factor of twelve by combining the four fold uniform Cartesian mask with a three fold pseudo-random variable density undersampling pattern.

To create the single channel data, we performed a principal component analysis (PCA) on the original multi-channel data and selected the most significant component to obtain a coil compressed data. To generate the undersampled Fourier measurements, we retrospectively undersampled the data using a 30% uniform random sampling mask.

1) data_2D_withoutmotion_complex.mat: This is the coil combined data (N1xN2xNt) obtained by doing a PCA on the original multi-channel data.
2) mask0_ur(0.3)forpm.mat: This mat file corresponds to the sampling mask whose size is equal to the size of the image series. The ones in each frame of the mask correspond to the sampled locations.
3) T2multicoildata.mat: This represents the multi-channel Fourier data (N1xN2xNtxNcoil).
4) truthforT2.mat: This mat file is obtained from a SENSE reconstruction using the multichannel data. The dimension is N1xN2xNt.
5) mask0var()forpm.mat: This mat file corresponds to the sampling mask, which is used for the multi-channel experiments.
6) mask_improved_dataset2.mat : This is another sampling mask which is used to mask the background and CSF regions in the image series; the maps are then estimated from the masked data.
