% This program plots the roots of
% 20 random degree 20 polynomials.
% (However, the plot here is renormalized 
% by dividing by the arctan of the norm 
% of the root, thus yielding the uniform 
% distribution on the unit disc.)
% The distribution of the coefficients
% is complex Gaussian with covariance
% matrix diagonal with entries square roots
% of binomial coefficients...
%    Sept 23, 2002, J. Maurice Rojas

hold off; newplot; 

bb=diag(binom(20,[0:20]));

axis equal; hold on; 

for j=1:20
 r=roots(bb*(randn(21,1)+i*randn(21,1)));
 plot(r.*atan(abs(r))./abs(r),'.'); 
end
