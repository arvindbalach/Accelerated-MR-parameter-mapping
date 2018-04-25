function x_recon = xsub_sc_direct(x_recon,Ahb,AtA,Q,res,lambda,filter_siz,p)
%%% Solves sub-problem 2 for the coil combined case: equation (19) from the paper.

f = @(x)(regul_direct(x,Q,filter_siz,res,p) + lambda*AtA(x));

rhs = lambda*Ahb;
[x_recon,flag1,res1,iter1,err1] = pcg(f,rhs(:),1e-8,25,[],[],x_recon(:));
x_recon = reshape(x_recon,res);

clear Ahb Q  

end