% todo: fix any bugs which may pop up. 
% questions and comments welcome! 
% J. Maurice Rojas, August 3, 2005 
% 
% usage: tot=amoeba(supp,coeffs,piclim,sig,res,nn,pp,fig,soup) 

function tot=amoeba(supp,coeffs,piclim,sig,res,nn,pp,fig,soup) 

% supp=2xm array of support vectors in Z^2
% coeffs=1xm array of coefficients
% piclim=[minx,maxx,miny,maxy], for picture... 
% sig= 0, 1, or 2, according as you want the positive part, the 
%              entire unfolded real part, or the folded real part of the amoeba...
% e.g., res=0; % draw just the complex part, no real part...
% e.g., res=400; % draw the real part too, fairly finely...
% nn = x-density of points (across half width of picture) 
% pp = phase-density used for each coordinate of base in fibered drawing...
% fig = figure number
% soup = string indicating support k-tuple...

tot=0; 

% # of steps for real part is proportional 
%  to largest of width/height...

eps=.0001; % tolerance for imag<eps => real...

figure(fig);
hold off; newplot; axis equal; axis xy; axis off; hold on; 

 % draw complex part...
 % fibered over x-axis...
 sot=fiber(supp,coeffs,piclim,nn,pp,1,[1 1],0,'b.'); tot=tot+sot;
 % fibered over y-axis...
 sot=fiber(supp,coeffs,piclim,nn,pp,2,[1 1],0,'b.'); tot=tot+sot;
    
axis equal; axis on; 

% draw the complex or ``positive'' amoeba, according as doreal is 0 or not,   
% fibered over the x_coord axis... furthermore, flip either the picture  
% or the signs, according as doreal is 1 or 2...
function tot=fiber(supp,coeffs,piclim,nn,pp,coord,flip,doreal,col); 

 m=size(supp,2); tot=0; 

 % sort the exponents for later evaluation along coord-axes...
 sm=max(supp(3-coord,:)); 
 s1=min(supp(3-coord,:)); 
 d=sm-s1; 

 for i=0:(nn-1) % index for norm of zi
  for p=0:(pp-1) % index for phase of zi 
   tot=tot+1;
   nzi=piclim(2*coord-1)+i*(piclim(2*coord)-piclim(2*coord-1))/(nn-1);
   pzi=p*2*pi/pp; 
   if (doreal~=0)&(flip(coord)==-1) 
    pzi=pzi+pi; 
   end; 
   zi=exp(nzi+sqrt(-1)*pzi);
   c=zeros(1,d+1);

   % looks like this simplified version works better! 
   % no plot mistakes so far! 8/31/05 
   for k=1:m 
    c(d+1+s1-supp(3-coord,k))=c(d+1+s1-supp(3-coord,k))+coeffs(k)*zi^supp(coord,k);
   end;
   r=roots(c)'; nzjs=log(abs(r)); imj=abs(imag(r)); 
   nzj=[]; t=0; ss=size(nzjs,2); 

   % dirty trick: if you're doing the real part, 
   %   then filter out the roots with big imag 
   %   parts and real parts with the wrong sign...
   %   to force this, make the norm huge...
   if doreal~=0
    imj=imj+10*(flip(3-coord)*real(r)<=0); 
    nzjs=nzjs+2*(piclim(6-2*coord)-piclim(5-2*coord))*(imj>=eps);
   end; 

   for k=1:ss 
    if (nzjs(k)<=piclim(6-2*coord))&(nzjs(k)>=piclim(5-2*coord))
     nzj=[nzj nzjs(k)]; t=t+1;
    end;
   end;
   nzi=nzi*ones(1,t); 
   if doreal~=2  
    if coord==1
     plot(flip(1)*nzi,flip(2)*nzj,col);
     % plot(flip(1)*(nzi-piclim(1)),flip(2)*(nzj-piclim(3)),col);
    else 
     plot(flip(1)*nzj,flip(2)*nzi,col);
     % plot(flip(1)*(nzj-piclim(1)),flip(2)*(nzi-piclim(3)),col);
    end;
   else 
    if coord==1
     plot(nzi,nzj,col);
     % plot(nzi-piclim(1),nzj-piclim(3)),col);
    else 
     plot(nzj,nzi,col);
     % plot(nzj-piclim(1),nzi-piclim(3)),col);
    end;
   end;  
  end;
 end; 
