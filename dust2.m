% This program plots the roots of
% 2 random degree 50 polynomials.
% The distribution of the coefficients
% is complex Gaussian with covariance
% matrix diagonal with entries square roots  
% of binomial coefficients... 
%    Sept 23, 2002, J. Maurice Rojas

hold off; newplot; 

bb=diag(binom(50,[0:50]));

axis equal; hold on; 

for j=1:2
 plot(roots(bb*(randn(51,1)+i*randn(51,1))),'.'); 
end


