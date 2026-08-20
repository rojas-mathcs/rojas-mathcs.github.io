% 3/17/09, copyright J. Maurice Rojas
% This function sets up the matrices needed to 
% draw a slice of Amoeba(Nabla_A), corresponding 
% either a coordinate subspace (determined by cell) 
% or the canonical slice (orthogonal to the row space of 
% Ahat). Usage: 
% function [a,cell,logvars,n,m,cano,nb,bmap,map,smap,shifts]=setslice(a,cell,logvars,quiet);  
% ...where 
%  a        = n x m support matrix
%  cell     = n-subset of [1,...,m] with det(a(:,cell)) odd...
%  logvars  = 1 x m array containing the desired logs of the variances
%              for the real gaussians will be using...
%   quiet   = flag for shutting up (1) or not (0)  

% and the output is
%  bmap  = extra transformation for walls...
%  map   = the main transformation sending the true amoeba to
%          (a) the sliced defined by cell
%               = -a_cell\a_noncell = (a(:,cell)^{-1) a(:,noncell)
%               (if cell does not contain 0),
%           OR
%           (b) the canonical slice (meaning the slice of the amoeba
%               orthogonal to the rows of Ahat
%
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

function [a,cell,logvars,n,m,cano,nb,bmap,map,smap,shifts]=setslice(a,cell,logvars,quiet);  

% set up the support and cell (i.e., the coeffs c1 and ck are
% homogenized down to 1, where k is the index in cell...)
n=size(a,1); m=size(a,2); 

% need size(logvars,2)=m
if size(logvars,2)~=m
  if quiet~=1 
   disp(sprintf('Your log-variance vector is of the wrong-dimension!\n Using zeroes for the logs to proceed...\n'));
  end; 
 logvars=zeros(1,m);
end;

% consider the canonical slice option...
if ismember(0,cell)
 cano=1; cell=0;
else
 cano=0;
end;

% sort exponents by last coordinate AND remember to reset cell AND
% remember to resort logvars!!!
% this sorts in the usual way in the univariate case and, in higher
% dimensions, keeps distinct polynomial supports separate if the
% cayley trick is used...
[y,ind]=sort(a(n,:));
[y,invy]=sort(ind); % get inverse of permutation defined by ind...
a=a(:,ind);
logvars=logvars(ind(:));
if cano==0
 cell=invy(cell(:));
end;

aa=a; % save the sorted, but UNshifted, version of the support A...

% force A to have first point = origin...
origif=(a(:,1)==zeros(n,1));
if sum(origif(:))<n
 a=a-a(:,1)*ones(1,m);
end;

% now work out the matrices that will effect the slicing map...
% a basis for the right nullspace of Ahat...
% note: using null without 'r' allows matlab to use
%       some orthgonality, which seems to yield better
%       numerical conditioning. also, in the projection
%       formulae that follow in the plotting subroutines,
%       I ASSUME AN ORTHONORMAL BASIS!!!
nb=null([ones(1,m);a]);
shifts=null([ones(1,m);a],'r'); % rational version for perusal...
if cano==0
 % to do the coordinate slice corresponding to cell, you need to know
 % -(A_C)^{-1} A_{C'}, which would just be -a(:,cell)\a(:,noncell)...
 noncell=setdiff([2:m],cell); % the complement of the odd cell in {2,...,m}
 pram=-a(:,cell)\a(:,noncell);
 pa=mod(a,2); % for determining signs later...
 spram=mod(-pa(:,cell)\pa(:,noncell),2); % NOT FAIL-SAFE HERE!:
                                    % YOU NEED THE CELL TO BE ODD!
                                    % i.e., odd determinant for a(:,cell),
                                    % assuming a(:,1)=origin...
 bmap=zeros(m,m-n-1); map=zeros(m,m-n-1); smap=zeros(m,m-n-1);
 bmap(noncell,:)=eye(m-n-1);
 map(noncell,:)=bmap(noncell,:);
 map(1,:)=-(ones(1,n)*pram+ones(1,m-n-1));
 map(cell,:)=pram;
 smap(noncell,:)=bmap(noncell,:);
 smap(1,:)=mod(-(ones(1,n)*spram+ones(1,m-n-1)),2);
 smap(cell,:)=spram;
else
 map=nb; % for rnearckthkplot.m 
 bmap=map; % for drawrays.m 
 smap=2*ones(m,m-n-1); % flipping quadrants might be iffy for canononical slice
end;

% revert A to its sorted but UNshifted version...
a=aa; 
