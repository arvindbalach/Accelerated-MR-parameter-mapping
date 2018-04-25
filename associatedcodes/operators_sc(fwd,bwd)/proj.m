function y = proj(x,ind)
%%% Computes AtA(x)
y = gpuArray.zeros(size(x));
y(ind) = x(ind);
y = y(:);
end