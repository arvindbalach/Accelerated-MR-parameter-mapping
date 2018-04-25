function result = AhAmc_gpu(x,overres,res,csm,Samp,ncoils,ind_full)
%%% Defines the AhA operator
res_ov = gpuArray.zeros(overres);
n1 = res(1);n2 = res(2); nf = res(3);
na = overres(1);nb = overres(2); nc = overres(3);
x = reshape(x,na,nb,nc);
x = x(ind_full);
x = ifft2(reshape(x,res));

result = gpuArray.zeros(n1*n2,nf);
Samp = reshape(Samp,n1*n2,nf);
for i=1:nf
    xi  = x(:,:,i);
    S_r = gpuArray.zeros(n1,n2);
    Si = Samp(:,i);
    ind = find(Si);
    clear Si;
    for j=1:ncoils
        temp2 = gpuArray.zeros(n1,n2);
                        tmp = fft2(csm(:,:,j).*xi);
%         tmp = fft2(csm(:,:,j,i).*xi);
        temp2(ind) = tmp(ind); %%%% measured samples onto zero filled matrix
                        S_r = S_r+(conj(csm(:,:,j)).*ifft2(temp2));
%         S_r = S_r+(conj(csm(:,:,j,i)).*ifft2(temp2));
        clear tmp temp2
    end
    S_r = fft2(S_r);
    result(:,i) = S_r(:);
    clear S_r
end
result = reshape(result,res);
res_ov(ind_full) = result;
result = res_ov(:);
end