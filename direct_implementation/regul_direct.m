function term = regul_direct(x,Q,filter_siz,res,p)
%%% regul_direct: Computes the gradient of the first term in equation (19)
n1 = res(1);
n2 = res(2);
n3 = res(3);
x =  gather(reshape(x,n1,n2,n3));
Q = gather(Q);

x_s = fftshift(fftshift(x,1),2);
[n1_s,n2_s,n3_s] = size(x_s);

%%% Forms Toeplitz matrix
Toep_x = im2colstep(real(x_s),filter_siz,[1,1,1]) + 1i*im2colstep(imag(x_s),filter_siz,[1,1,1]);
clear x;
Toep_x = rot90(Toep_x,-1);

temp= Toep_x*(Q*Q');

%%% Forms the Adjoint
tempr = rot90(temp);
Toep_adj = col2imstep(real(tempr),[n1_s,n2_s,n3_s],filter_siz) + 1i*col2imstep(imag(tempr),[n1_s,n2_s,n3_s],filter_siz);
Toep_adj = ifftshift(ifftshift(Toep_adj,1),2);

term = (2/p)*gpuArray(Toep_adj(:));

end
