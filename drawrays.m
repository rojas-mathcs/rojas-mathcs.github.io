% questions and comments welcome!
% J. Maurice Rojas; 2/8/09  
% This code sorts and draws the rays of the amoeba of 
% an A-discriminant for given A in Z^n of cardinality n+3 (using 
% some formulae I derived from the Horn-Kapranov uniformization). In 
% particular, this code draws the chamber cones by allowing adjacent 
% rays to intersect at a common origin.  
% usage: drawrays(nn,rays,cs,ineq,m,piclim,fold,tol);
% where...
%  nn      = REDUCTION of m x m matrix whose ith row N_i is the normal to 
%            the ith A-discriminant wall (and the corresponding shift just 
%            happens to be log(max(1,|N_i|) ) 
%            (so nn here should be an m x (m-n-1) matrix...) 
%  cs      = m x 1 matrix, consisting of the constant terms for the
%            equations ax+by=c describing the ith wall...
%  ineq    = REDUCTION of m x m matrix whose ith row I_i decides which 
%            half-plane defines the ith A-discriminant wall W_i via W_i 
%            being { x | N_i(x-s_i)=0 AND I_i.x>>0} 
%            (so ineq here should be an m x (m-n-1) matrix...) 
%  piclim  = [minx,maxx,miny,maxy], for picture
%  fold    = 1 iff you want the rays to be drawn for the UNfolded amoeba.  
%  tol     = tolerance for what will be declared to be zero...
function drawrays(nn,rays,cs,ineq,m,piclim,fold,tol); 

% set up the angles with the positive x-ray to sort...
for i=1:m 
 % note: i'm assuming that nn has NO zero row. 
 %       under this assumption, we need NOT check 
 %       whether (abs(rays(i,1))<tol)&(abs(rays(i,2))<tol) ... 
 if abs(rays(i,1))<tol % in principle, i could just check if this integer is zero...
  angle(i)=sign(rays(i,2))*pi/2; 
 elseif abs(rays(i,2))<tol % in principle, i could just check if this integer is zero...
  angle(i)=pi*(1-sign(rays(i,1)))/2; % 0 (resp. pi) if on pos (resp. neg) ray... 
 else 
  xangle=atan(rays(i,2)/rays(i,1)); 
  angle(i)=xangle-sign(xangle)*(pi*(1-sign(rays(i,1)))/2); 
 end; 
end; 
% find order of rays by angle...
[y,ind]=sort(angle); 

disp(sprintf('These are the oriented walls!')); 
rays 

% now sort everybody before proceeding...
nn=nn(ind(:),:); cs=cs(ind(:)); ineq=ineq(ind(:),:); rays=rays(ind(:),:); 

p=zeros(m,2); q=zeros(m,2); 
for i=1:m
 % now determine what kind of intersection to draw...
 % check intersection with left, right, bottom, and top...
 left=0; right=0; bottom=0; top=0; 
 % check left/right...
 if abs(nn(i,2))>tol 
  y=(cs(i)-nn(i,1)*piclim(1))/nn(i,2); 
  if (y>=piclim(3))&(y<=piclim(4))
   left=1; 
  end; 
  y=(cs(i)-nn(i,1)*piclim(2))/nn(i,2); 
  if (y>=piclim(3))&(y<=piclim(4))
   right=1; 
  end;
 end;  
 % check down/up...
 if abs(nn(i,1))>tol
  x=(cs(i)-nn(i,2)*piclim(3))/nn(i,1);
  if (x>=piclim(1))&(x<=piclim(2))
   bottom=1;
  end;
  x=(cs(i)-nn(i,2)*piclim(4))/nn(i,1);
  if (x>=piclim(1))&(x<=piclim(2))
   top=1;
  end;
 end;

 done=0; 
 if left==1 
  done=done+1; 
  p(i,:)=[piclim(1) (cs(i)-nn(i,1)*piclim(1))/nn(i,2)]; 
 end; 
 if right==1
  done=done+1; 
  if done==2 
   q(i,:)=[piclim(2) (cs(i)-nn(i,1)*piclim(2))/nn(i,2)]; 
  else 
   p(i,:)=[piclim(2) (cs(i)-nn(i,1)*piclim(2))/nn(i,2)]; 
  end;  
 end; 
 if (bottom==1)&(done<2) 
  done=done+1;
  if done==2 
   q(i,:)=[(cs(i)-nn(i,2)*piclim(3))/nn(i,1) piclim(3)];
  else
   p(i,:)=[(cs(i)-nn(i,2)*piclim(3))/nn(i,1) piclim(3)];
  end;
 end;
 if (top==1)&(done<2)
  done=done+1;
  if done==2 
   q(i,:)=[(cs(i)-nn(i,2)*piclim(4))/nn(i,1) piclim(4)];
  else
   p(i,:)=[(cs(i)-nn(i,2)*piclim(4))/nn(i,1) piclim(4)];
  end;
 end;

 % now figure out how to plot the rays!: check which 
 % extreme is farther in the correct direction, then connect to 
 % the intersection of 2 adjacent rays...
 if (p(i,:)-q(i,:))*ineq(i,:)'>0 
  q(i,:)=p(i,:); 
 end; 
 % intersection with next ray... 
 pf(i,:)=[nn(i,:); nn(1+mod(i,m),:)]\[cs(i); cs(1+mod(i,m))];
 % intersection with previous ray... 
 pb(i,:)=[nn(i,:); nn(1+mod(i-2,m),:)]\[cs(i); cs(1+mod(i-2,m))];

 % plot ray segments from bounding box intersection to the forward and 
 % backward cone vertices...
 plot([pf(i,1) q(i,1)]-piclim(1),[pf(i,2) q(i,2)]-piclim(3),'r-'); 
 plot([pb(i,1) q(i,1)]-piclim(1),[pb(i,2) q(i,2)]-piclim(3),'r-'); 
 % note: i need both intersections, for otherwise i have missing 
 %       ray pieces (observed 2009/01/07) 
end; 
