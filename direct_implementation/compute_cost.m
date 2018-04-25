function cost = compute_cost(p,s,eps,x,b,A,lambda0,cost)

if p == 0
    schatten = 0.5*sum(log(s-eps));
else
    schatten = (1/p)*sum((s-eps).^(p/2));
end
% 
% figure(100);plot(gather(0.5*log10(abs(s))));drawnow;
% if p == 0
%     schatten = 0.5*sum(log(abs(s)-eps));
% else
%     schatten = (1/p)*sum(abs((s)-eps).^(p/2));
% end
diff = A(x)-b;
thiscost = lambda0*schatten + 0.5*norm(diff(:)).^2; %objective function

clear diff x b;

cost = [cost,thiscost];
figure(1); plot(gather(abs(cost)),'-k*'); title('cost'); drawnow;