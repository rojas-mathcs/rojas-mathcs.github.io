% questions and comments welcome!
% J. Maurice Rojas; 2/8/09  
% This code illustrates the discriminant chambers and 
% chamber cones (on log paper) for the family of polynomials 
% supported on any near-circuit (subset A of Z^n with cardinality n+3),  
% by using my implementation of horn-kapranov uniformization and 
% some additional tricks (rnearckthkplot.m).  
%
% More precisely, this code draws the real 1-dimensional part of the 
% underlying A-discriminant variety, and the chamber cone walls, 
% on log paper. This program also takes the trouble to sort the 
% walls so that the cones are known. 
% 
% usage: 
% [a cell n m nb bmap map smap shifts]=chambers(a,cell,logvars,piclim,dots,ray,res,fig) 
% where 
%  a       = n x m support matrix
%  cell    = n-subset of [1,...,m] with det(a(:,cell)) odd...
%  logvars = 1 x m array containing the desired logs of the variances
%  piclim  = [minx,maxx,miny,maxy], for picture 
%  dots    = 1 or 0 according as you want some random polynomials drawn or not 
%  ray     = 1 or 0 according as you want the chamber walls drawn or not 
%  res     = plot resolution = number of points drawn per curve branch...
%  fig     = figure number where drawing will appear
% and the output is 
%  nb    = basis for right null space of ahat
%  bmap  = extra transformation for walls... 
%  map   = the main transformation sending the true amoeba to
%          (a) the sliced defined by cell
%               = -a_cell\a_noncell = (a(:,cell)^{-1) a(:,noncell) 
%               (if cell does not contain 0), 
%           OR 
%           (b) the canonical slice (meaning the slice of the amoeba 
%               orthogonal to the rows of Ahat 
%  smap  = mod 2 version of map that use for exponentiation to find 
%          the sign of the orthant you should draw in         
%  shifts  = rational basis for right nullspace for ahat
%          = null([ones(m);a],'r')
% ...so the true discriminant variety is nothing more than
%  {u.* [t^a(:,:)] | ahat*u=0 , t in Csn}
%  or
%  {u.* [t^a(:,:)] | u in colspace(shift), t in Csn}
% and the reduced discriminant amoeba drawn is nothing more than
% (a) 
% {log( u_noncell/u_1.*u_cell/u_1^pram ) | u in colspace(shift) }
% or
% { log(u_cell/u_1)*pram + log(u_noncell/u_1) | u in colspace(shift) }
% OR 
% (b) 
% {log( u*map ) | u in colspace(shift) } 
function [a cell n m nb bmap map smap shifts]=chambers(a,cell,logvars,piclim,dots,ray,res,fig) 

% get the matrices you'll need to plot your slice of Amoeba(Nabla_A), 
% and sort a, cell, and logvars along the way...
[a cell logvars n m cano nb bmap map smap shifts]=setslice(a,cell,logvars,0); 

% set up for plotting, if the dimension is low enough!
if m==n+3
 figure(fig); hold off; newplot; axis xy; hold on;
 lowerleftcorner=[piclim(1) piclim(3)];
 % set up title, according to whether you're peppering with random
 % polynomials or not...
 polystring=setdiscplottitle(a,n,m,cell,cano,dots);
 if dots==1 
  numpts=100; % a default value for the number of random polynomials to use...
  title(sprintf(polystring,numpts,logvars));
  % pepper with random polynomials, if requested...
  pepperhorn(m,map,numpts,logvars,lowerleftcorner);
 else
  numpts=[]; logvars=[];
  title(sprintf(polystring));
 end;
end; 

% force A to have first point = origin...
origif=(a(:,1)==zeros(n,1));
if sum(origif(:))<n
 disp(sprintf('Warning!: Things go better when the first point of A is O, so I''m shifting your support...'));
 a=a-a(:,1)*ones(1,m);
end;

disp(sprintf('Your sorted and shifted support is...'));
disp(a); 

ahat=[ones(1,m); a]; 
[nn,cs,walls,ineq]=findwalls(ahat,n,m,nb,bmap); 

% if the dimension is low enough, go ahead and plot! 
if m==n+3

 % don't flip any orthants, set the resolution of the curve drawing, 
 % fold everything on to the positive quadrant, set number of random 
 % tetranomials to pepper with, set variances, 
 % set which figure you'll use, and set tolerance for what zero is... 
 flip=[1 1]; fold=2; tol=1e-10;  
 
 rnearckthkplot(a,nb,map,smap,flip,piclim,res,fold,tol); 

 % finally, if desired, compute and draw the amoeba rays
 if ray==1
  drawrays(nn,walls,cs,ineq,m,piclim,fold,tol); 
 end; 

 % throw in the projection of the origin in R^m (which happens to 
 % just be the origin in R^{m-n-1}), just for fun...
 plot(0-lowerleftcorner(1),0-lowerleftcorner(2),'gx'); 

 % readjust the window, just in case... 
 axis([0 piclim(2)-piclim(1) 0 piclim(4)-piclim(3)]);

 % the plot looks prettier without the axes, and they're
 % a bit intricate to describe concisely and correctly anyway...
 axis equal; axis off; 

end; 
