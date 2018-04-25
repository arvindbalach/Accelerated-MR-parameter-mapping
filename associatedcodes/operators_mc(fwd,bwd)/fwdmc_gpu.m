function A_mc_gpu = fwdmc_gpu(x,csm,ind,n1,n2,nf,ncoils)
%%%% The Forward operator is defined;Here x is the Fourier data

temp = gpuArray.zeros(n1,n2,ncoils,nf);
x = reshape(x,n1,n2,nf);
factor=1;
I = ifft2(x);
for i = 1:nf
    xi = I(:,:,i);
    for j = 1:ncoils
        temp(:,:,j,i) = factor*fft2(csm(:,:,j).*xi);
    end
end
A_mc_gpu = temp(ind);
clear temp x

end