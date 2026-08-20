% questions and comments welcome!
% J. Maurice Rojas; 5/13/09  
% This is an example of what the real part of a certain A-discriminant 
% looks like, along with an illustration of its chamber cones. 
% In particular, this example is mentioned in Examples 2.7, 2.10, 
% 2.12, and 3.2 of [BHPR10]. 

% set support and cell...
a=[0 404 405 808]; cell=[3]; 

% set the bounding box and logs of the variances for the real centered gaussians... 
piclim=[-1 2 -3 -.2]; 
logvars=[0 8 8 0]; 
% note: log(808 choose 404) is about 556.4895401 and
%       log(808 choose 405) is about 556.4870679...

% no dots now, but please draw the chamber walls...
dots=0; rays=1; 

% high resolution...
res=500; 

% do the picture on figure 1...
chambers(a,cell,logvars,piclim,dots,rays,res,fig); 

% find point near nabla_A...
ll=[piclim(1) piclim(3)]; 
plot(0-ll(1),0-ll(2),'gx'); 

title('');  
