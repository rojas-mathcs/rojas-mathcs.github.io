% questions and comments welcome!
% J. Maurice Rojas; 3/18/09  
% This code finds the walls (hyperplanes) defined by the membranes of 
% Amoeba(Nabla_A).  
% underlying polynomial f, thus finding the canonical lifting 
% for f when f lies in an unbounded chamber. 
% The underlying technique is to do redundancy checks of constraints 
% (via LP) and then check (via LP again) which walls are unbounded. 
% If you have at least n walls that are unbounded in the correct 
% direction, then you MUST have a unique chamber cone containing your point, 
% and you get the cone vertex (aka the proper shift for ArchNewt!) as 
% output.   
% usage:  [nn,cs,walls,ineq]=findwalls(ahat,n,m,nb,bmap);  
%  where 
%   ahat     = (n+1) x m LIFTED support matrix
%   nb       = right null space of ahat 
%   bmap     = extra linear transformation for the walls... 
%   map      = the main linear transformation sending the true amoeba to
%               (a) the amoeba of the reduced discriminant defined by cell
%                   = -a_cell\a_noncell = -a(:,cell)^{-1} a(:,noncell)
%                   (if cell does not contain 0)
%               OR
%               (b) the canonical slice of the amoeba, orthogonal to
%                   the row space of Ahat
function [nn,cs,walls,ineq]=findwalls(ahat,n,m,nb,bmap)  

% % ith row = ith inequality vector
% % so this is a stacking of
% % ahat*ahat' [alpha;beta] = ahat * e_i
% % (and you can also make the matrix integral,
% % if you like, by throwing in afactor of
% % det(ahat*ahat')...)
ineq=ahat'*((ahat*ahat')\ahat)-eye(m); 
% ineq=-eye(m); % even easier: just consider how the    
              % parametrization blows up near the poles!!!! 
ineq=ineq*bmap; % no project the inequality vector... 
                                             
% now, to check membership, let's prepare the reduced normals and 
% reduced inequality vectors... (pretty much as in drawrays...) 
% first, let's do the normals and shifts (similar to how they were 
% done in drawrays).  

if m==n+2 
 % this is the best way for n+3...
 nn=nb*[0 -1; 1 0]*nb'; % new slick formula discovered on 3/20/09!
                        % WAY faster than old det method!: over 279
                        % times faster when n=500 (!)
 ss=log(abs(nn)); % ith row = ith shift
 for i=1:m
  ss(i,i)=0;
 end;
 % in particular, since the entries of nb may be NON-integral,
 % you can no longer just say log(max(1,ss))...
else 
 % this is the easy way to generalize...
 % NOTE: here, i'm assuming that all consecutive 
 %       (m-n-2)-tuples of rows of B are linearly independent! 
 nn=zeros(m); ss=zeros(m); k=m-n-2; 
 for i=1:m
    
  pole=null(nb(1+mod([0:(k-1)]+i-1,m),:)); 

  ss(i,:)=log(abs((nb*pole)')); 
  ss(i,:); ss(i,1+mod([0:k-1]+i-1,m))=zeros(1,k); 
    
  stack=[zeros(k,i-1) ones(k,1) zeros(k,m-i)]; 
  for j=1:k 
   stack(j,1+mod([1:(k-1)]+i+j-2,m))=ones(1,k-1);
  end;

  nn(i,:)=null([stack;ahat])';

 end; 
end; 

% here are the ``constants'' on the righthand side of
% ax+by=c...
for i=1:m
  cs(i)=nn(i,:)*ss(i,:)'; 
end;

% now we reduce the normals...
nn=nn*bmap; 

% and figure out the directed rays, when the dimension is low enough !!!! 
if m==n+3 
 walls=[nn(:,2) -nn(:,1)]; 
 for i=1:m
  if walls(i,:)*ineq(i,:)'<0 
   walls(i,:)=-walls(i,:); 
  end; 
 end; 
else 
 walls=[]; % future code for m!=n+3 here...
end; 
