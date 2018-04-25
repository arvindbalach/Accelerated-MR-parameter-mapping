%%
% Code for the publication:
%%% Arvind Balachandrasekaran, Vincent Magnotta and Mathews Jacob, "Recovery of damped exponentials
%%% using structured low rank matrix completion" Vol 36,2087-2098, Oct 2017.
%%% Authors: Arvind Balachandrasekaran, Vincent Magnotta and Mathews Jacob

%%% Date:  Nov22, 2017

%%% Main file for reconstructing coil-combined kspace data from
%%% undersampled Fourier measurements.
%%% We solve the following optimization problem:min_x ||Ax-b||^2 + lambda0||Toep(x)||_p; x is the Fourier data to be
%%% recovered and ||X|_p refers to the Schatten p norm of X, where 0<=p<=1
%%% We introduce Fourier domain approximations to solve eq (17) and (18) in the paper.
%% Add paths.
clear all;close all;clc
set(0,'DefaultFigureWindowStyle','docked');
addpath('../associatedcodes/direct-liftandunlift-codes/')
addpath('../associatedcodes/operators_sc(fwd,bwd)/')
addpath('../associatedcodes/getkspaceindices')
addpath('../data/');
%% 
%%%% Load the data
load data_2D_withoutmotion_complex.mat
load mask0_ur(0.3)forpm.mat
X0 = data_complex(:,:,13:end); %%% corresponds to T2 weighted images.
mask  = mask0;
X0 = gpuArray(X0);
m_abs= max(abs(X0(:)));
X0 = X0/m_abs;
[n1,n2,nf] = size(X0);
m = (1/sqrt(n1*n2))*fft2(X0);
%% Set filter sizes and define other associated variables 
res=[n1 n2 nf];
f1 = 7;f2 =7;ft=11;
filter_siz = [f1 f2 ft];
filter_siz2 = [2*f1 2*f2 ft] - [1,1,0];
overres = res +0*[2*f1 2*f2 0];
na = overres(1);nb = overres(2);nc = overres(3);
%%
%%% Get k space indices for filter and the kspace data
[ind_full,ind_filter] = get_kspaceindices_proposed(overres,res,filter_siz);
k1 = get_kspace_inds([overres(1),overres(2)]);
ind_filter2d = get_lowpass_inds(k1,[filter_siz2(1) filter_siz2(2)]);
%%
ind_samples = find(mask0~=0);
[A,At] = defAAt(ind_samples,res);
b = A(m);
mask_pad = gpuArray.zeros(overres);
mask_pad(ind_full) = mask0;
ind_samples_pad = find(mask_pad~=0);
AtA = @(x)proj(x,ind_samples_pad);
Atb = At(b);
Atb_pad = gpuArray.zeros(overres);
Atb_pad(ind_full) = Atb; %%% they are the same when res = overres.
x = Atb;
x_pad = gpuArray.zeros(overres);
x_pad(ind_full) = x;
%% 
%%% Define the optimization parameters.
p=0.6;%% p for schatten p-norm
q = 1-(p/2);
lambda0 =(8e-3)*p/2;  %regularization parameter (high value enforces low-rank constraint)
lambda = 1/lambda0;
eps = initial_epsilon(x_pad,0,filter_siz,filter_siz2,ind_filter2d,overres);
eta=1.4;
epsmin = 1e-9;   %minimum possible epsilon value
%%
%%% Define other variables such as cost etc.
cost = [];
SNR = [];
Xr = (sqrt(n1*n2))*ifft2(x);%%% Ground truth image
[mu,s]= estimate_mu(x_pad,eps,filter_siz,filter_siz2,ind_filter,ind_filter2d,overres,q);
cost_t0 = compute_cost(p,s,eps,x,b,A,lambda0,[]);%%% Initial cost
SNR_t0 = 20*log10(norm(X0(:))/norm(Xr(:)-X0(:)));%%% Initial SNR
ti = [];
tic;
maxiter = 45; %%% Iterations can be reduced.
%%
%%%% Solve equation (17) and (18) in an alternating fashion, stop when 
%%%% convergence or max iterations reached
for i=1:maxiter
    i
    xtemp =x_pad;
    [mu,s]= estimate_mu(xtemp,eps,filter_siz,filter_siz2,ind_filter,ind_filter2d,overres,q);
    ChC = pre_computeChC(mu,filter_siz,res);
    x_pad = xsub_hybrid_sc(x_pad,Atb_pad,AtA,ChC,overres,lambda,p);
    x = reshape(x_pad(ind_full),res);
    cost = compute_cost(p,s,eps,x,b,A,lambda0,cost);
    eps = max(eps/eta,epsmin);
    Xr = (sqrt(n1*n2))*ifft2(x);
    SNR_iter = 20*log10(norm(X0(:))/norm(Xr(:)-X0(:)));
    SNR = [SNR SNR_iter];
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

TE = [10:10:120]';
[T2g,~] = t2_MapsCalc(x1.*repmat(maskforT2,[1,1,12]),TE);
[T2p,~] = t2_MapsCalc(x2.*repmat(maskforT2,[1,1,12]),TE);

T2g_cr = T2g;
T2p_cr = T2p;

SNRp=20*log10(norm(T2g_cr(:))/norm(T2p_cr(:)-T2g_cr(:)))
SNR(end)
% keyboard