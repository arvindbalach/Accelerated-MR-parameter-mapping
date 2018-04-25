function ChC = pre_computeChC(mu,filter_siz,res)
%%% Computes G corresponding to equation (26) in the paper.
%%%% For some details, please refer to the paragraph describing the formation of G in the
%%%% paper.
fc = filter_siz(3);
num = size(mu,4);
T = res(3);
k = T-fc+1;
ChC = gpuArray.zeros(res(1)*T);
na = res(1);nb = res(2);
for i=1:num
    M = gpuArray.zeros(res(1)*T);
    M_sum = gpuArray.zeros(res(1)*T);
    %      mu1 = gpuArray(mu(:,:,:,i));
    mu1 = mu(:,:,:,i);
    
    %%% forms mu1_res and mu2_res without using for loops.
    mu1_temp = permute(conj(mu1),[1,3,2]);
    mu1_res = reshape(mu1_temp,[],size(mu1,2),1);
    mu2_res = reshape(mu1,size(mu1,2),[],1);
    %%
    %%% Forms mu1_res and mu2_res using for loops.
    %    tic;
%         l=1;
%         for j=1:fc
%             mu1_res(l:l+res(1)-1, 1:res(2)) = conj(mu1(:,:,j));
%             mu2_res(1:res(1),l:l+res(2)-1) = mu1(:,:,j);
%             l=l+res(1);
%         end
%         clear l j;
    
    mu1_rep = repmat(mu1_res,[1,fc]);
    mu2_rep = repmat(mu2_res,[fc,1]);
    mu_out = mu1_rep.*mu2_rep;
    
    M(1:na*fc,1:nb*fc) = mu_out;
    clear mu1 mu1_temp mu1_res mu2_res mu1_rep mu2_rep mu_out
    for count = 1:k
        M_sum =M_sum + M;
        M = circshift(M,[res(1),res(2)]);
    end
    clear count M;
    ChC = ChC + M_sum;
    clear M_sum
end
end