function epsinit= initial_epsilon(x,epsilon,filter_siz)
%%% epsinit initializes the epsilon value for the IRLS minimization.
%%% It is defined as 0.01*max(singular value^2 of the Toeplitz matrix
%%% formed from the zero filled kspace data.
%%
%%%% The im2colstep uses cpu variables.
x = gather(x);
x = fftshift(fftshift(x,1),2);
R = im2colstep(real(x),filter_siz,[1,1,1]) + 1i*im2colstep(imag(x),filter_siz,[1,1,1]);
clear x;

R = rot90(R,-1); %%% converts to Toeplitz
gram_R = R'*R; 
[~,S] = eig(gpuArray(gram_R)+epsilon*gpuArray.eye(size(gram_R)));
clear R gram_R;
s = diag(S); %note: s = sing. values squared.
epsinit = 0.01*max(s(:)); %%% make sure this is real valued, if not use abs(s(:))
end