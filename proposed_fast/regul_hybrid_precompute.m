function res = regul_hybrid_precompute(x,ChC,overres,p)
%%%% This functions computes the first term of equation (26) from the
%%%% paper. This is the first term of the gradient of the modified cost function in
%%%% (25) in the paper.
vec = @(z)(z(:));
const1 = sqrt(overres(1)*overres(2));
I =  const1*ifft2(reshape(x,overres));
clear x
T = overres(3);
k=1;
res = gpuArray.zeros(overres);
for i=1:T
%     temp = reshape(ChC(k:k+128-1,1:end),128,128,T);
    temp = reshape(ChC(k:k+overres(1)-1,1:end),overres(1),overres(2),T);
    res(:,:,i) = sum((temp.*I),3);
    clear temp
%     k =k+128;
    k =k+overres(1);
end
clear I
res = (2/p)*(1/const1)*vec(fft2(res));
end