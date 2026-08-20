% Copyright 2017, J. Maurice Rojas
% 
% This function converts numbers to symbols, 
% as requested in Problem 1b on HW#1. 
% In particular, if you save this file as nc.m in 
% a directory, run Matlab in the same directory,  
% and then say 
%  nc([24 21 28 28 31])
% you should get 
%  ans =
% 
%  hello 
% This code is merely meant to serve as an example 
% of a suitable solution. 

function st=nc(nu) 

l=length(nu);  
for i=1:l 
 if nu(i)==2 % space 
  st(i)=' ';
 elseif nu(i)==3  
  st(i)='''';
 elseif nu(i)==4  
  st(i)=',';
 elseif nu(i)==5  
  st(i)='-';
 elseif nu(i)==6 % . 
  st(i)='.'; 
 elseif (nu(i)>6)&(nu(i)<17) 
  st(i)=char(nu(i)+41);  
 else 
  st(i)=char(nu(i)+80); 
 end; 
end; 
