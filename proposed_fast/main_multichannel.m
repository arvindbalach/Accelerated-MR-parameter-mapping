% Code for the publication:
%%% Arvind Balachandrasekaran, Vincent Magnotta and Mathews Jacob, "Recovery of damped exponentials
%%% using structured low rank matrix completion" Vol 36,2087-2098, Oct 2017.
%%% Authors: Arvind Balachandrasekaran, Vincent Magnotta and Mathews Jacob

%%% Date:  Nov22, 2017

%%% Main file for reconstructing multi-channel kspace data from
%%% undersampled Fourier measurements.
%%% We solve the following optimization problem:min_x ||Ax-b||^2 + lambda0||Toep(x)||_p; x is the Fourier data to be
%%% recovered and ||X|_p refers to the Schatten p norm of X, where 0<=p<=1
%%% We introduce Fourier domain approximations to solve eq (17) and (18).

%% Add paths.
clear all;close all;clc
set(0,'DefaultFigureWindowStyle','docked');
addpath('../associatedcodes/direct-liftandunlift-codes/')
addpath('../associatedcodes/operators_mc(fwd,bwd)/')
addpath('../associatedcodes/getkspaceindices')
addpath('../data/');
%%
%%%% Load the data
load T2multicoildata.mat
load truthforT2.mat;
load csmforT2.mat
load mask0(var3*4).mat
data = gpuArray(kdata);
[n1,n2,ncoils,nf] = size(data);
mask  = mask0;
csm = gpuArray(csm);
m_abs= max(abs(X0(:)));
X0 = X0/m_abs;
data = data/m_abs;
%% Set model order
res = [n1 n2 nf];
f1 = 27;f2 = 27;ft=3;
filter_siz = [f1 f2 ft];
filter_siz2 = [2*f1 2*f2 ft] - [1,1,0];
overres = res + 0*[2*f1 2*f2 0];
na = overres(1);nb = overres(2);nc = overres(3);
%%
%%% Get k space indices for filter and the kspace data
[ind_full,ind_filter] = get_kspaceindices_proposed(overres,res,filter_siz);
k1 = get_kspace_inds(overres);
ind_filter2d = get_lowpass_inds(k1,[filter_siz2(1) filter_siz2(2)]);
%%
Samp = gpuArray(mask);
masktemp= repmat(Samp,[1,1,1,ncoils]);
for j=1:ncoils
    mask_allcoils(:,:,j,:)  = masktemp(:,:,:,j);
end
clear masktemp;
ind = find(mask_allcoils(:));
clear mask_allcoils
%%
conjcsm = conj(csm);
A_mc_gpu = @(x)fwdmc_gpu(x,csm,ind,n1,n2,nf,ncoils);
Ah_mc_gpu = @(z)bwdmc_gpu(z,conjcsm,ind,n1,n2,nf,ncoils);
b = data(ind);
clear data;
Ahb = Ah_mc_gpu(b);
x = Ahb;
%%
%%% Define the optimization parameters.
p=0.7;%%% p for schatten p-norm
q = 1-(p/2);
lambda0 =5*p/2;  %regularization parameter (high value enforces low-rank constraint)
lambda = 1/lambda0;
eps = initial_epsilon(x,0,filter_siz,filter_siz2,ind_filter2d,overres);
eta=1.4;
epsmin = 1e-9;   %minimum possible epsilon value
%%
%%% Define other variables such as cost etc.
cost = [];
SNR = [];
ti = [];
tic;
maxiter = 20;
%%
%%%% Solve equation (17) and (18) in an alternating fashion, stop when
%%%% convergence or max iterations reached
for i=1:maxiter
    i
    [mu,s]=estimate_mu(x,eps,filter_siz,filter_siz2,ind_filter,ind_filter2d,overres,q);
    ChC = pre_computeChC(mu,filter_siz,res);
    x=xsub_hybrid_mc(x,Ahb,ChC,csm,Samp,overres,res,lambda,ncoils,ind_full,p);
    cost = compute_cost(p,s,eps,x,b,A_mc_gpu,lambda0,cost);
    clear mu;
    Xr = (sqrt(n1*n2))*ifft2(x);
    max(abs(Xr(:)))   
    eps = max(eps/eta,epsmin);
    SNR_iter = 20*log10(norm(X0(:))/norm(Xr(:)-X0(:))); 
    SNR = [SNR,SNR_iter];
    figure(2); plot(gather(SNR)); title('SNR'); drawnow;
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
keyboard
