function epsinit= initial_epsilon(x_pad,epsilon,filter_siz,filter_siz2,ind_filter2d,overres)
%%% eps_proposed initializes the epsilon value for the IRLS minimization.
%%% It is defined as 0.01*max(eigen value)  of the gram matrix
%%% formed from the zero filled kspace data.
%%
%%%% Gram matrix formed according to equation (27), (28) and (29)
fa = filter_siz(1);fb = filter_siz(2);fc =filter_siz(3);
na = overres(1);nb = overres(2);nc = overres(3);
xtemp = gather(x_pad);
clear x_pad;
const1=1;
count=1;
for k=fc:-1:1
    temp = xtemp(:,:,k:k+nc-fc);
    ftemp = conj(const1*ifft2(temp));
    clear temp;
    count1 = 1;
    for j=fc:-1:1
        temp1 = xtemp(:,:,j:j+nc-fc);
        gr1 = sum((1/const1)*fft2(const1*ifft2(temp1).*ftemp),3);
        gr_sh = fftshift(reshape(gr1(ind_filter2d),filter_siz2(1),filter_siz2(2)));
        R((count-1)*fa*fb + 1:count*fa*fb,(count1-1)*fa*fb + 1:count1*fa*fb) = rot90(im2colstep(real(gr_sh),[fa,fb],[1 1])+1i*im2colstep(imag(gr_sh),[fa,fb],[1 1]),-1);
        clear gr_sh gr1;
        count1=count1+1;
    end
    count = count+1;
end

R = (na*nb)*R;
[~,S] = eig(gpuArray(R)+epsilon*gpuArray.eye(size(R)));
clear R;
s = diag(S); %note: s = sing. values squared.
epsinit = gather(abs(0.01*max(s(:))));
end