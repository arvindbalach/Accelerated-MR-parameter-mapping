function Ah_mc_gpu = bwdmc_gpu(z,conjcsm,ind,n1,n2,nf,ncoils)
%%%% The backward operator is defined;

temp = gpuArray.zeros(n1,n2,ncoils,nf);
temp(ind) = z;
Ah_mc_gpu = gpuArray.zeros(n1,n2,ncoils,nf);
factor = 1;
for i = 1:nf
    for j = 1:ncoils
        Ah_mc_gpu(:,:,j,i) = (factor*ifft2(temp(:,:,j,i))) .*conjcsm(:,:,j);
    end
end
Ah_mc_gpu = squeeze(sum(Ah_mc_gpu,3));
Ah_mc_gpu = fft2(Ah_mc_gpu);
end