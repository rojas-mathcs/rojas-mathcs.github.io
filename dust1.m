% This program plots the roots of 
% 30 random degree 20 polynomials.  
% The distribution of the coefficients 
% is complex Gaussian with covariance 
% matrix the identity matrix...
%    Sept 23, 2002, J. Maurice Rojas
hold off; newplot; 

axis equal; hold on; 

for j=1:30
 plot(roots(randn(21,1)+i*randn(21,1)),'.'); 
end

