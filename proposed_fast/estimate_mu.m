function [mu,s]= estimate_mu(xtemp,eps,filter_siz,filter_siz2,ind_filter,ind_filter2d,overres,q)
%%% Estimates mu which is used to form the D matrix; this matrix is then
%%% used in the second sub-problem.
%%
%%% Gram matrix formed according to equation (27), (28) and (29)
fa = filter_siz(1);fb = filter_siz(2);fc =filter_siz(3);
na = overres(1);nb = overres(2);nc = overres(3);
xtemp = gather(xtemp);
const1 = sqrt(na*nb);
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
        count1=count1+1;
        clear gr1 gr_sh;
    end
    count = count+1;
end
R = (na*nb)*R; %%%% just to match the scale with direct implementation.
%%
%%% eigen decomposition of the gram matrix R
[Utemp,Stemp] = eig(R+eps*eye(size(R)));
clear R;
stemp = diag(Stemp); %note: s = sing. values squared.
clear Stemp;
[s,ind] = sort(stemp,'ascend');
U = Utemp(:,ind);
clear S Utemp stemp

sqrt_alpha = s.^(-q/2);
Q= U*diag(sqrt_alpha);
clear U;
%%
%%% Estimates the polynomials mu corresponding to each column of Q
%%%% For very large filters, if there is no gpu memory convert to cpu
%%%% variables and then estimate the polynomials.
const2 = sqrt(na*nb);
mu = gpuArray.zeros(na,nb,fc,length(abs(s)));
for i = 1:length(abs(s))
    h = gpuArray.zeros(overres);
    h(ind_filter) = ifftshift(reshape(Q(:,i),filter_siz));
    %%% At this point, the frames of the filter are of the form :
    %%% h3 h4 h5 h1 h2 (if the number of frames is 5)
    temp = const2*ifft2(h);
    clear h;
    %%% The following if else code does the following: The filter frames are as follows:
    %%% they are [ h5 h4 h3 h2 h1]
    if mod(fc,2)
        mu_temp = temp(:,:,[fliplr((1:(fc-1)/2 +1)),fliplr((end-(fc-1)/2 +1 :end))]);
    else
        mu_temp = temp(:,:,[fliplr((1:fc/2)),fliplr((end-(fc/2) +1 :end))]);
    end
    clear temp;
    mu(:,:,1:fc,i) = mu_temp;
    clear mu_temp;
end
end