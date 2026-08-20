% Copyright 2017, J. Maurice Rojas
% 
% This function converts symbols to numbers, 
% as requested in Problem 1a on HW#1. 
% In particular, if you save this file as cn.m in 
% a directory, run Matlab in the same directory,  
% and then say 
%  cn('hello')
% you should get 
%  ans =
% 
%    24    21    28    28    31 
% This code is merely meant to serve as an example 
% of a suitable solution. 

function nu=cn(s) 

l=length(s); nu=double(s); 
for i=1:l 
 if nu(i)==32 % space 
  nu(i)=2;
 elseif nu(i)==39 % ' 
  nu(i)=3;
 elseif nu(i)==44 % , 
  nu(i)=4;
 elseif nu(i)==45 % - 
  nu(i)=5;
 elseif nu(i)==46 % . 
  nu(i)=6; 
 elseif (nu(i)>46)&(nu(i)<58) 
  nu(i)=nu(i)-41;  
 else 
  nu(i)=nu(i)-80; 
 end; 
end; 
