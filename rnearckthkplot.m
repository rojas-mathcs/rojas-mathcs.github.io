%  questions and comments welcome!
% J. Maurice Rojas; 3/7/09  
% This code plots the 1-dimensional real part of an A-discriminant 
% variety (using the Horn-Kapranov uniformization along with some 
% extra implementation tricks), on log paper, in the special case that 
% A is in Z^n and #A=n+3. 
% usage: rnearckthkplot(a,nb,map,smap,flip,piclim,res,fold,tol);
%  where...                                
%  a          = n x m support matrix (assumed to have a(:,1)=origin!!!)
%                                     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%  nb         = ORTHONORMAL basis for right null space of Ahat := 
%               [ones(1,size(a,2)); a];  
%  map        = the main linear transformation sending the true amoeba to 
%               (a) the amoeba of the reduced discriminant defined by cell 
%                   = -a_cell\a_noncell = -a(:,cell)^{-1} a(:,noncell) 
%                   (if cell does not contain 0) 
%               OR 
%               (b) the canonical slice of the amoeba, orthogonal to 
%                   the row space of Ahat 
%  smap       = cousin of map that helps keep track of orthant signs 
%               by allowing one to plot the UNfolded amoeba slice via 
%               !!!!!  
%  flip       = possible flip of signs, for simplification/visualization 
%               purposes...
%  piclim     = [minx,maxx,miny,maxy], for picture
%  res        = density of points along real direction 
%  fold       = plot pos part, unfolded real part, or folded real part,
%               according as fold=0,1,2
%       tol   = tolerance for what will be declared to be zero...
%
% ...so that the true discriminant variety is nothing more than 
%  {u.* [t^a(:,:)] | ahat*u=0 , t in Csn}  
%  or 
%  {u.* [t^a(:,:)] | u in colspace(nb), t in Csn}  
% so that the reduced discriminant amoeba drawn is nothing more than 
% {log( u_noncell/u_1.*u_cell/u_1^pram ) | u' in colspace(nb) } 
% or
% { log(u_cell/u_1)*pram + log(u_noncell/u_1) | u' in colspace(nb) } 

% i don't think i need to worry about division by zero any more. 
% so let's shut if off for now...
warning off MATLAB:divideByZero  % doesn't work on my older matlab versions...

function rnearckthkplot(a,nb,map,smap,flip,piclim,res,fold,tol);

m=size(a,2); n=size(a,1);  

% where are the corresponding poles, in semicircular coordinates? 
poles=sort(-atan(nb(:,1)./nb(:,2))'); np=m; % number of poles...
poles=[-.5*pi poles .5*pi]; np=size(poles,2); 

% put in the coordinate cross, thickened to indicate ambiguity, 
% if you are unfolding...
if fold==1
 thick=.1;  
 % yellow and green ambiguity cross...
 fill([piclim(1)-piclim(2) piclim(2)-piclim(1) piclim(2)-piclim(1) piclim(1)-piclim(2)],thick*[piclim(3)-piclim(4) piclim(3)-piclim(4) piclim(4)-piclim(3) piclim(4)-piclim(3)],'y');  
 plot([piclim(1)-piclim(2) piclim(2)-piclim(1) piclim(2)-piclim(1) piclim(1)-piclim(2) piclim(1)-piclim(2)],thick*[piclim(3)-piclim(4) piclim(3)-piclim(4) piclim(4)-piclim(3) piclim(4)-piclim(3) piclim(3)-piclim(4)],'y-');  
 fill(thick*[piclim(1)-piclim(2) piclim(2)-piclim(1) piclim(2)-piclim(1) piclim(1)-piclim(2)],[piclim(3)-piclim(4) piclim(3)-piclim(4) piclim(4)-piclim(3) piclim(4)-piclim(3)],'y');  
 plot(thick*[piclim(1)-piclim(2) piclim(2)-piclim(1) piclim(2)-piclim(1) piclim(1)-piclim(2) piclim(1)-piclim(2)],[piclim(3)-piclim(4) piclim(3)-piclim(4) piclim(4)-piclim(3) piclim(4)-piclim(3) piclim(3)-piclim(4)],'y-');  
 fill([piclim(1)-piclim(2) piclim(2)-piclim(1) piclim(2)-piclim(1) piclim(1)-piclim(2)],.1*thick*[piclim(3)-piclim(4) piclim(3)-piclim(4) piclim(4)-piclim(3) piclim(4)-piclim(3)],'g');
 plot([piclim(1)-piclim(2) piclim(2)-piclim(1) piclim(2)-piclim(1) piclim(1)-piclim(2) piclim(1)-piclim(2)],.1*thick*[piclim(3)-piclim(4) piclim(3)-piclim(4) piclim(4)-piclim(3) piclim(4)-piclim(3) piclim(3)-piclim(4)],'g-');
 fill(.1*thick*[piclim(1)-piclim(2) piclim(2)-piclim(1) piclim(2)-piclim(1) piclim(1)-piclim(2)],[piclim(3)-piclim(4) piclim(3)-piclim(4) piclim(4)-piclim(3) piclim(4)-piclim(3)],'g');
 plot(.1*thick*[piclim(1)-piclim(2) piclim(2)-piclim(1) piclim(2)-piclim(1) piclim(1)-piclim(2) piclim(1)-piclim(2)],[piclim(3)-piclim(4) piclim(3)-piclim(4) piclim(4)-piclim(3) piclim(4)-piclim(3) piclim(3)-piclim(4)],'g-');
end; 
 
% now let's do the real part!
for k=1:(np-1)    

 if poles(k)<poles(k+1)  % don't bother unless the interval is non-empty...

  % figure out which quadrant the image of your segment is going to lie in...
  mid=(poles(k)+poles(k+1))/2; 
  signell=sign(nb*[cos(mid);sin(mid)])';  

  if fold~=2 
   amosign=[prod(signell.^(smap(:,1)')) prod(signell.^(smap(:,2)'))] 
   plotsign=flip.*amosign; 
  else 
   plotsign=[1 1]; 
  end; 

  % adjust the signs, according to whether 
  % you're plotting just the positive part, 
  % the unfolded real part, or the folded real part...
  if ((plotsign(1)==1)&(plotsign(2)==1))|(fold==1) 
   for j=1:(res-2)  
    t1=poles(k)+j*(poles(k+1)-poles(k))/res; t2=t1+(poles(k+1)-poles(k))/res;  
    % we'll need the log absolute values, and the signs, of the 
    % various linear forms... 
    % ...where nullbasis=nb as above... 
    logabsell1=log(abs(nb*[cos(t1);sin(t1)]))'; logabsell2=log(abs(nb*[cos(t2);sin(t2)]))'; 

     amoabs1=logabsell1*map; amoabs2=logabsell2*map;  

    if ((amoabs1(1)>=piclim(1))&(amoabs1(1)<=piclim(2))&(amoabs1(2)>=piclim(3))&(amoabs1(2)<=piclim(4)))|((amoabs2(1)>=piclim(1))&(amoabs2(1)<=piclim(2))&(amoabs2(2)>=piclim(3))&(amoabs2(2)<=piclim(4)))  
     plot(plotsign(1)*([amoabs1(1) amoabs2(1)]-piclim(1)),plotsign(2)*([amoabs1(2) amoabs2(2)]-piclim(3)),'b-'); 
    end; % if you're in the plot window...
   end; % for j...
  end; % if the signs are correct...
  % adjust the window, according to whether
  % you're plotting just the positive part,
  % the unfolded real part, or the folded real part...
  if (fold==0)|(fold==2)
   axis([0 piclim(2)-piclim(1) 0 piclim(4)-piclim(3)]); 
  else % fold had better be 1 at this point...
   axis(2*piclim); 
  end;
 end; % if subinterval is OK...
end; % for k over interval indices

% adjust the window, according to whether
% you're plotting just the positive part,
% the unfolded real part, or the folded real part...
if (fold==0)|(fold==2)
 axis([0 piclim(2)-piclim(1) 0 piclim(4)-piclim(3)]);
else % fold had better be 1 at this point...
 axis([piclim(1)-piclim(2) piclim(2)-piclim(1) piclim(3)-piclim(4) piclim(4)-piclim(3)]);
end;

