
# Structured matrix completion algorithm for accelerated Magnetic Resonance (MR) parameter mapping
# PUBLICATION
A. Balachandrasekaran, V. Magnotta and M. Jacob. "Recovery of damped exponentials using structured low rank matrix completion." *IEEE Transactions on Medical Imaging*, 36(10):2087-2098, 2017. 
# Codes

The main files inside the proposed_fast and direct_implementation folders solve the following optimization problem

![](https://latex.codecogs.com/gif.latex?%24%5Cmin_%7B%5Cwidehat%20%5Crho%7D%20%5C%7CT%28%5Cwidehat%20%5Crho%29%5C%7C_p%20&plus;%20%28%5Cmu/2%29%5C%7CA%28%5Cwidehat%20%5Crho%29-b%5C%7C%5E2_2%24)

with and without introducing approximations respectively. ![](https://latex.codecogs.com/gif.latex?%24T%28%5Cwidehat%20%5Crho%29%24) is the Toeplitz matrix formed from the Fourier samples ![](https://latex.codecogs.com/gif.latex?%24%5Cwidehat%20%5Crho%24), ![](https://latex.codecogs.com/gif.latex?%24%5C%7CX%5C%7C_p%24) is the Schatten p-norm of ![](https://latex.codecogs.com/gif.latex?%24X%24), ![](https://latex.codecogs.com/gif.latex?%24A%24) is the Fourier undersampling operator and ![](https://latex.codecogs.com/gif.latex?%24b%24) represents undersampled Fourier measurements.
