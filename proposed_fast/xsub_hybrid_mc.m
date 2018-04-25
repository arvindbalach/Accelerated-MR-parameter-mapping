function x_pad = xsub_hybrid_mc(x_pad,Ahb_pad,ChC,csm,mask,overres,res,lambda,ncoils,ind_full,p)
%%% Solves sub-problem 2 for the multichannel case: equation (19) from the paper.
f = @(x)(regul_hybrid_precompute(x,ChC,overres,p)+ lambda*AhAmc_gpu(x,overres,res,csm,mask,ncoils,ind_full));
rhs = lambda*Ahb_pad;
[x_pad,flag1,res1,iter1,err1] = pcg(f,rhs(:),1e-8,10,[],[],x_pad(:));
clear Ahb_pad mu mask_pad csm 
x_pad = reshape(x_pad,overres(1),overres(2),overres(3));
end