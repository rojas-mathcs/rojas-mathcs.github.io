function p=swaprow(A,i,j); 

p=A; 
s=p(i,:); 
p(i,:)=p(j,:); 
p(j,:)=s; 
