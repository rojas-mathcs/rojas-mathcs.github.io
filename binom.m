% This program computes the binomial 
% coefficient ``n choose i'' by computing 
% Pascal's triangle row by row...
%   Sept. 23, 2002, J. Maurice Rojas
function b=binom(n,i);

a=zeros(1,n+1);
aa=zeros(1,n+1);

a(1)=1;
a(2)=1;
aa(1)=1;

if i==0
 b=1;
elseif n==0
 b=0;
elseif n==1
 b=1;
elseif or((i>n),(i<0))
 b=0;
else
 for j=2:n
  for k=1:j
   aa(k+1)=a(k+1)+a(k);
  end;
  aa(j+1)=1;
  a=aa;
 end;
 b=a(i+1);
end

