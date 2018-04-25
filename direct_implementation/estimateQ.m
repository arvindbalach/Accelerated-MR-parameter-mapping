function [Q,s]= estimateQ(x,eps,filter_siz,q)
%%% Solves sub-problem 1: estimateQ estimates the matrix Q, which is required for solving eqn (18)
%%%% Refer section weight update in the paper.
x = gather(x);
x = fftshift(fftshift(x,1),2);
R = im2colstep(real(x),filter_siz,[1,1,1]) + 1i*im2colstep(imag(x),filter_siz,[1,1,1]);
clear x
R =rot90(R,-1);
gram_R = R'*R;
[U,S] = eig(gram_R+eps*eye(size(gram_R)));
s = diag(S); %note: s = sing. values squared.
clear S;
sqrt_alpha = s.^(-q/2);
Q= U*diag(sqrt_alpha);
end