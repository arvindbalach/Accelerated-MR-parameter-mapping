%%
%%% Code by Arvind Balachandrasekaran
%%% Date:  Nov22, 2017

%%% Main file for reconstructing coil-combined kspace data from
%%% undersampled Fourier measurements.
%%% We solve the following optimization problem:min_x ||Ax-b||^2 + lambda0||Toep(x)||_p; x is the Fourier data to be
%%% recovered and ||X|_p refers to the Schatten p norm of X, where 0<=p<=1
%%% There are no approximations introduced here and we implement equations
%%% (17) and (18) from the paper* directly.
%%% For details refer paper*: A.Balachandrasekaran et al., "Recovery of damped exponentials
%%% using structured low rank matrix completion".
%% Add paths.
clear all;close all;clc
set(0,'DefaultFigureWindowStyle','docked');
addpath('../associatedcodes/direct-liftandunlift-codes/')
addpath('../associatedcodes/operators_sc(fwd,bwd)/')
addpath('../data/');
%%
%%%% Load the data
load data_2D_withoutmotion_complex.mat
load mask0_ur(0.3)forpm.mat
X0 = data_complex(:,:,13:end);
mask  = mask0;
X0 = gpuArray(X0);
m_abs= max(abs(X0(:)));
X0 = X0/m_abs;
[n1,n2,nf] = size(X0);
m = (1/sqrt(n1*n2))*fft2(X0);
%% Set filter sizes and define other associated variables 
res = [128,128,12];
f1 = 7;f2 = 7;ft=11; %%% This is an optimization parameter
filter_siz = [f1 f2 ft];
%% 
ind_samples = find(mask0~=0);
[A,At] = defAAt(ind_samples,res);
b = A(m);
AtA = @(x)proj(x,ind_samples);
Atb = At(b);
x = Atb;
%%
%%% Define the optimization parameters.
p=0.6; %% p for schatten p-norm
q = 1-(p/2);
lambda0 = (5e-7)*p/2;  %regularization parameter (high value enforces low-rank constraint)
lambda = 1/lambda0;
eps= initial_epsilon(x,0,filter_siz); %%% epsilon initialization:
eta=1.4;
epsmin = 1e-9;   %minimum possible epsilon value
%%
%%% Define other variables such as cost etc.
cost = [];
term1=[]; 
term2=[];
SNR = [];

Xr = (sqrt(n1*n2))*ifft2(x); %%% Ground truth image
[~,s]= estimateQ(x,eps,filter_siz,q);
cost_t0 = compute_cost(p,s,eps,x,b,A,lambda0,[]);%%% Initial cost
SNR_t0 = 20*log10(norm(X0(:))/norm(Xr(:)-X0(:)));%%% Initial SNR

ti = [];
tic;
maxiter = 45;
%%
%%%% Solve equation (17) and (18) in an alternating fashion, stop when 
%%%% convergence or max iterations reached
for i=1:maxiter
    i
    xtemp =x;
    [Q,s]= estimateQ(xtemp,eps,filter_siz,q);
    x = xsub_sc_direct(x,Atb,AtA,Q,res,lambda,filter_siz,p);
    cost = compute_cost(p,s,eps,x,b,A,lambda0,cost);
    eps = max(eps/eta,epsmin);
    Xr = (sqrt(n1*n2))*ifft2(x);
    SNR_iter = 20*log10(norm(X0(:))/norm(Xr(:)-X0(:)));
    SNR = [SNR,SNR_iter];
    figure(2),plot(SNR);title('SNR');drawnow;
    ti = [ti toc];
    if i>1
        if abs(cost(end)-cost(end-1))/abs(cost(end))<=1e-6
            break;
        end
    end
end
Xp = gather(Xr);
X0 = gather(X0);
%%
%%%% Computing the T2 maps.

load mask_improved_dataset2.mat
maskforT2 = imrotate(maskforT2,-90);

x1 = imrotate(X0,-90);
x2 = imrotate(Xp,-90);

TE = [10:10:120]'; %%% Echo times
[T2g,~] = t2_MapsCalc(x1.*repmat(maskforT2,[1,1,12]),TE);
[T2p,~] = t2_MapsCalc(x2.*repmat(maskforT2,[1,1,12]),TE);

T2g_cr = T2g;
T2p_cr = T2p;

SNRp=20*log10(norm(T2g_cr(:))/norm(T2p_cr(:)-T2g_cr(:)))
SNR(end)