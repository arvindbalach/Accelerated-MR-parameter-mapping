function x_pad = xsub_hybrid_sc(x_pad,Ahb_pad,AtA,ChC,overres,lambda,p)
%%% Solves sub-problem 2 for the coil combined case: equation (19) from the paper.
f = @(x)(regul_hybrid_precompute(x,ChC,overres,p) + lambda*AtA(x));
rhs = lambda*Ahb_pad;
[x_pad,flag1,res1,iter1,err1] = pcg(f,rhs(:),1e-8,25,[],[],x_pad(:));
x_pad = reshape(x_pad,overres);
clear Ahb_pad mu csm_pad mask_pad 
end