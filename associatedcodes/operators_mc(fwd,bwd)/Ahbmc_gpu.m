function Ahb = Ahbmc_gpu(b,conjcsm,ind,n1,n2,nf,ncoils)

%%% Another function to compuite Ahb; We could just use bwdmc_gpu function.
%%% The result from this function and the bwdmc_gpu is the same.
temp = gpuArray.zeros(n1,n2,ncoils,nf);
res = gpuArray.zeros(n1*n2,nf);
temp(ind)  = b;
factor=1;
for i=1:nf
    tmp = gpuArray.zeros(n1,n2,ncoils);
    for j=1:ncoils
        tmp(:,:,j) = conjcsm(:,:,j).*(factor*ifft2(temp(:,:,j,i)));
    end
    tmp = fft2(tmp);
    res(:,i)= reshape(sum(tmp,3),n1*n2,1);
end
Ahb = res(:);
end

