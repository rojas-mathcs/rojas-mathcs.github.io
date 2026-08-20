function p=pivot(A,i,j); 

s=size(A,1); 
p=eye(s); 
p(:,i)=-A(:,j)/A(i,j); 

%%%%%%%%%% To make it pivoting from diag down...
for k=1:(i-1)
 p(k,i)=0;
end; 
%%%%%%%%%%%%%%

p(i,i)=1/A(i,j); 
