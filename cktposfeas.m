% Copyright, J. Maurice Rojas; 2/11/09  
% (Any use of this code within other software must acknowledge the 
% author. Use within commercial software is permitted only with 
% the author's consent.) 
% This code decides positive feasibility for honest n-variate 
% (n+2)-nomials (checking whether the inputs are honest beforehand), 
% via Algorithm 3.2 from [Bihan-Rojas-Stella, 2009].  
%
% Warning: This implemenation does NOT use gcd-free bases or 
%          high-brow linear forms in logs estimates: it merely 
%          checks linear forms in logs via Matlab's arithmetic. 
%          Nevertheless, this implementation has worked well on 
%          numerous examples.  Comments and questions are welcome: 
%          rojas@math.tamu.edu 
%
% usage: 
% cktposfeas(a,coeffs) 
% where 
%  a       = n x m support matrix
%  coeffs  = 1 x m coefficient array 
%     tol  = small positive number for zero threshold... 
%  ...and a defines a circuit (degenerate allowed).  (In essence, 
%  this means that a is the union of a point with the vertex set 
%  of a simplex.) The output is a (hopefully true) declaration as 
%  to whether f(x) := c_1 x^{a_1} + c_2 x^{a_2} + ... + c_m x^{a_m} has a 
%  positive root or not.  The complexity is polynomial in n+log(deg(f)), 
%  unlike older methods that were polynomial in deg(f)^n.   
function cktposfeas(a,coeffs,tol) 

% set up the support 
n=size(a,1); m=size(a,2); 

% throw out supports that are too large...
if m>n+2 
 disp(sprintf('Your support has too many points!')); 
 return; 
end; 

% keep track of coefficient signs...
posc=[]; pc=0; 
zerc=[]; zc=0; 
negc=[]; ngc=0; 
for j=1:m 
 if coeffs(j)>0
  posc=[posc j]; pc=pc+1;
 elseif abs(coeffs(j))<tol
  zerc=[zerc j]; zc=zc+1;
 else
  negc=[negc j]; ngc=ngc+1;
 end;
end; 

% throw out inputs with zero coefficients...
if zerc>0 
 disp(sprintf('Zero coefficients not allowed! Please re-enter your polynomial in the correct format.')); 
end; 

% finish early if there are no sign alternations...
if (pc==0)|(ngc==0) 
 disp(sprintf('Your polynomial has NO positive roots since it''s coefficients are all of the same sign. Duh!'));  
 return; 
end; 

% unique difference of sign made negative...
if pc==1 
 coeffs=-coeffs; 
 ngc=pc; 
 negc=posc; 
end; 

% compute unique circuit relation (or see that you do not have a circuit)...
ahat=[ones(1,m); a];
b=null(ahat,'r'); 
if size(b,2)>1 
 disp(sprintf('Your support is similar to an (n''+k)-set in R^{n''} for k>2, and is tthus NOT a circuit...')); 
 return; 
end; 

% now let's sort and use my fancy algorithm...
% sort exponents AND coefficients using the b_i as keys. in particular, 
% we group by positive b_i, negative b_i, then vanishing b_i...
pos=[]; p=0; 
zer=[]; z=0; 
neg=[]; ng=0; 
for j=1:m
 if b(j)>0 
  pos=[pos j]; p=p+1;  
 elseif abs(b(j))<tol
  zer=[zer j]; z=z+1; 
 else 
  neg=[neg j]; ng=ng+1; 
 end; 
end; 

% unique difference of sign made negative...
if p==1
 b=-b;
 swap=ng; ng=p; p=swap; 
 swap=neg; neg=pos; pos=swap; 
end;

ind=[pos neg zer]; 
[y,inv]=sort(ind);     % find how indices will change...
b=b(ind(:));           % shuffle circuit relation coords...
neg=inv(neg(:));     % shuffle => change indices of negative circuit rel coords...
coeffs=coeffs(ind(:)); % same shuffle for coeffs...
negc=inv(negc(:));     % shuffle => change indices of negative coeffs...

if (ng==1)&(ngc==1)&(neg==negc)  
 jp=neg; % j' from my algorithm in [BRS09]... 
 % log(left-hand side of relations of Steps (b) & (c), Algor 3.2, pg 7, [BRS09]) 
 le=log([-b(jp) coeffs(1:(jp-1))])*[-b(jp);b(1:(jp-1))]; 
 % log(right-hand side of relations of Steps (b) & (c), Algor 3.2, pg 7, [BRS09]) 
 ri=log([-coeffs(jp) b(1:(jp-1))'])*[-b(jp);b(1:(jp-1))]; 

 if (jp<m)&(abs(le-ri)<tol)  
  % given the accruing sign conditions, equality would imply that 
  % the zero set intersects a coordinate hyperplane (after a monomial change 
  % of variables) in exactly a single point and does not meet the positive 
  % orthant...   
  disp(sprintf('Your polynomial has NO positive roots, but has a unique root on an orbit of positive codimension.')); 
  return;
 elseif le>ri   
  disp(sprintf('Your polynomial has NO positive roots, and stably so.')); 
  return; 
 end; 
end; 

%
disp(sprintf('Your polynomial has at least one positive root!')); 
