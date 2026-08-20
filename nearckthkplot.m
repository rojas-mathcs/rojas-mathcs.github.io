% questions and comments welcome!
% J. Maurice Rojas; 4/5/2010   
% This code plots the 1-dimensional real part of an A-discriminant 
% variety (using the Horn-Kapranov uniformization along with some 
% extra implementation tricks), on log paper, in the special case that 
% A is in Z^n and #A=n+3. 
% There are also options to plot the full amoeba and/or flipping 
% across orthants to recover sign information for the real part. 
% 
% usage: [pram shifts]=nearckthkplot(a,cell,flip,piclim,res,cxres,fold,cx,phasecolor,clr,tol);
%  where...                                
%  a          = n x m support matrix (assumed to have a(:,1)=origin!!!)
%  cell       = n-subset of [1,...,m] with det(a(:,cell)) odd...
%  flip       = possible flip of signs, for simplification/visualization 
%               purposes...
%  piclim     = [minx,maxx,miny,maxy], for picture
%  res        = density of points along real direction 
%  cxres      = density of points along phase direction 
%  fold       = plot pos part, unfolded real part, or folded real part,
%               according as fold=0,1,2
%  cx         = plot the traditional complex amoeba too!
%  phasecolor = vary (or not) the interior color according to the phase
%               of the pre-image under the Horn-Kapranov uniformization,
%               according as 1 or 0
%  col        = rgb vector for the color of the filled part of the amoeba...
%         tol = tolerance for what will be declared to be zero...
%
% ...and the output is 
%   pram   = the main transformation sending the true amoeba to 
%            the amoeba of the reduced discriminant defined by cell 
%          = -a_cell\a_noncell = -a(:,cell)^{-1} a(:,noncell) 
%  shifts  = rational basis for right nullspace for ahat 
%          = null([ones(m);a],'r') 
% ...so that the true discriminant variety is nothing more than 
%  {u.* [t^a(:,:)] | ahat*u=0 , t in Csn}  
%  or 
%  {u.* [t^a(:,:)] | u in colspace(shift), t in Csn}  
% so that the reduced discriminant amoeba drawn is nothing more than 
% {log( u_noncell/u_1.*u_cell/u_1^pram ) | u in colspace(shift) } 
% or
% { log(u_cell/u_1)*pram + log(u_noncell/u_1) | u in colspace(shift) } 
%
% examples: 
% (1) a=[0 1 2 3]; cell=2 ; flip=[1 1]; piclim=[-9 9 -9 9]; res=200; cxres=10; 
%     fold=2; cx=0; phasecolor=0; tol=1e-10;    
%  ...gives you the real part of the cubic discriminant.
% (2) a=[6 0 0 0 3 1; 0 3 1 6 0 0 ; 0 0 0 1 1 1 ]; cell=[3 4 6]; flip=[1 1]; 
%     piclim=[-6 6 -6 6]; res=800; cxres=20; fold=2; cx=1; phasecolor=1; 
%     tol=1e-10;    
%  ...gives you a plot of the amoeba of the discriminant of the 
%     rusek-shih family, augmented with some phase information.

% i don't think i need to worry about division by zero any more. 
% so let's shut if off for now...
warning off MATLAB:divideByZero  % doesn't work on my older matlab versions...

function [pram shifts]=nearckthkplot(a,cell,flip,piclim,res,cxres,fold,cx,phasecolor,clr,tol);

m=size(a,2); pa=mod(a,2); n=size(a,1);  
noncell=setdiff([2:m],cell); 
% the complement of the odd cell in {2,...,m}...

% force A to have first point = origin...
origif=(a(:,1)==zeros(n,1));
if sum(origif(:))<n
 disp(sprintf('Warning!: First point of a must be O. So I''m shifting your support...'));
 a=a-a(:,1)*ones(1,m);
end;

% to do the reduction, you need to know -(A_C)^{-1} A_{C'}, 
% which would just be -a(:,cell)\a(:,noncell)... 
% pram=-[2 3]; % the cubic works! 
% spram=[0 1]; % the cubic works!
% pram=[-33/35 12/35; -12/35 -2/35; 12/35 -33/35];  % rusek-shih works!
% spram=[1 0; 0 0 ; 0 1 ]; % tells us which signs matter...
       % entries are derived from parity of numerators 
       % of entries of pram. better: just do above formula MOD 2!
pram=-a(:,cell)\a(:,noncell); 
spram=mod(-pa(:,cell)\pa(:,noncell),2); % NOT FAIL-SAFE HERE!
                                        % YOU NEED AN ``ODD CELL''
                                        % i.e., odd determinant for a(:,cell)

% a basis for the right nullspace of A...
% note: using null with 'r' allows matlab to use 
%       some orthgonality, which seems to yield better 
%       numerical conditioning...
nb=null([ones(1,m);a]); 
shifts=null([ones(1,m);a],'r'); % output a rational version for perusal...

% set up for the plot...
figure(1); hold off; subplot(1,1,1); newplot; hold on;

% where are the corresponding poles, in semicircular coordinates? 
poles=sort(-atan(nb(:,1)./nb(:,2))'); np=m; % number of poles...
poles=[-.5*pi poles .5*pi]; np=size(poles,2); 

if cx==1 % if you want the traditional amoeba...
 for k=1:(np-1)  
   for j=1:(cxres-1) % here, I can hit the poles...
                                                                                
    t1=poles(k)+j*(poles(k+1)-poles(k))/cxres; t2=t1+(poles(k+1)-poles(k))/cxres;
    c1=cos(t1); s1=sin(t1);
    c2=cos(t2); s2=sin(t1);
                                                                                
    for an=1:(cxres-2) % leave a little gap so you don't hit the poles...

     % NOTE: Multiplying by the phase is NOT working well...   
     %       ...but going additively works REALLY well!... 
     %       in particular, adding an imag to the cos works OK, 
     %        although there are some weird holes...
     %       on the other, multiplying just the cos by the phase 
     %       seems to work best...
     ph1=.5*pi*(1+(sign(an-(cxres/2))*abs(an-(cxres/2))^.75/((cxres/2)^.75)));
     ph2=.5*pi*(1+(sign(an+1-(cxres/2))*abs(an+1-(cxres/2))^.75/((cxres/2)^.75))); 
% interesting: uniformizing via the above ugly ad hoc method 
%              actually improves the plot quite a bit!
%              the key difference is that, deeper in the interior, 
%              i don't have an unnecessarily high number of tiles... 11/2/06
%              ...even better: it improves my cubic plot too! 
%    ph1=pi*an/cxres; 
%     ph2=pi*(an+1)/cxres; 
     % (0,pi] appears to be the right range, and i can 
     %  probably prove it rigourously by modelling P^1_C 
     %  as a half-open hemi-sphere... 
     if phasecolor==1 
      clr=[sqrt(abs(sin(ph1/2))) abs(cos(ph1/2)) 0]; % check!!!
     end;  
        
     logabsella=log(abs(nb*[cos(t1)*exp(i*ph1);sin(t1)]))';
      logabsellb=log(abs(nb*[cos(t1)*exp(i*ph2);sin(t1)]))';
      logabsellc=log(abs(nb*[cos(t2)*exp(i*ph2);sin(t2)]))';
      logabselld=log(abs(nb*[cos(t2)*exp(i*ph1);sin(t2)]))';
         
     logabscellreda=logabsella(cell)-logabsella(1);
      logabscellredb=logabsellb(cell)-logabsellb(1);
      logabscellredc=logabsellc(cell)-logabsellc(1);
      logabscellredd=logabselld(cell)-logabselld(1);
           
     logabsnoncellreda=logabsella(noncell)-logabsella(1);
      logabsnoncellredb=logabsellb(noncell)-logabsellb(1);
      logabsnoncellredc=logabsellc(noncell)-logabsellc(1);
      logabsnoncellredd=logabselld(noncell)-logabselld(1);
 
     va=logabsnoncellreda+logabscellreda*pram;
      vb=logabsnoncellredb+logabscellredb*pram;
      vc=logabsnoncellredc+logabscellredc*pram;
      vd=logabsnoncellredd+logabscellredd*pram;
 
     if ((va(1)>=piclim(1))&(va(1)<=piclim(2))&(va(2)>=piclim(3))&(va(2)<=piclim(4)))|((vb(1)>=piclim(1))&(vb(1)<=piclim(2))&(vb(2)>=piclim(3))&(vb(2)<=piclim(4)))|((vc(1)>=piclim(1))&(vc(1)<=piclim(2))&(vc(2)>=piclim(3))&(vc(2)<=piclim(4)))|((vd(1)>=piclim(1))&(vd(1)<=piclim(2))&(vd(2)>=piclim(3))&(vd(2)<=piclim(4)))
      fill(([va(1) vb(1) vc(1) vd(1)]-piclim(1)),([va(2) vb(2) vc(2) vd(2)]-piclim(3)),clr); 
      % plot(([va(1) vb(1) vc(1) vd(1) va(1)]-piclim(1)),([va(2) vb(2) vc(2) vd(2) va(2)]-piclim(3)),'Color',clr); % change back to y- as necessary!!!
      plot(([va(1) vb(1) vc(1) vd(1) va(1)]-piclim(1)),([va(2) vb(2) vc(2) vd(2) va(2)]-piclim(3)),[clr '-']); % change back to y- as necessary!!!
  
      if fold==1 % if you needed unfolding...
       fill((-([va(1) vb(1) vc(1) vd(1)]-piclim(1))),([va(2) vb(2) vc(2) vd(2)]-piclim(3)),clr); 
       % plot((-([va(1) vb(1) vc(1) vd(1) va(1)]-piclim(1))),([va(2) vb(2) vc(2) vd(2) va(2)]-piclim(3)),'Color',clr); % change back to y- as necessary!!!
       plot((-([va(1) vb(1) vc(1) vd(1) va(1)]-piclim(1))),([va(2) vb(2) vc(2) vd(2) va(2)]-piclim(3)),[clr '-']); % change back to y- as necessary!!!
       fill((-([va(1) vb(1) vc(1) vd(1)]-piclim(1))),(-([va(2) vb(2) vc(2) vd(2)]-piclim(3))),clr);  
       % plot((-([va(1) vb(1) vc(1) vd(1) va(1)]-piclim(1))),(-([va(2) vb(2) vc(2) vd(2) va(2)]-piclim(3))),'Color',clr); % change back to y- as necessary!!!
       plot((-([va(1) vb(1) vc(1) vd(1) va(1)]-piclim(1))),(-([va(2) vb(2) vc(2) vd(2) va(2)]-piclim(3))),[clr '-']); % change back to y- as necessary!!!
       fill(([va(1) vb(1) vc(1) vd(1)]-piclim(1)),(-([va(2) vb(2) vc(2) vd(2)]-piclim(3))),clr); 
       % plot(([va(1) vb(1) vc(1) vd(1) va(1)]-piclim(1)),(-([va(2) vb(2) vc(2) vd(2) va(2)]-piclim(3))),'Color',clr); % change back to y- as necessary!!! 
       plot(([va(1) vb(1) vc(1) vd(1) va(1)]-piclim(1)),(-([va(2) vb(2) vc(2) vd(2) va(2)]-piclim(3))),[clr '-']); % change back to y- as necessary!!! 
      end; % if you needed unfolding...
     end; % if the square was in the plot limit...
    end; % for loop for angles...
   end; % for j...
 end; % for every subinterval... 
end; % if you want the traditional amoeba...

% put in the coordinate cross, thickened to indicate ambiguity,
% if you are unfolding...
if fold==1
 thick=.1;
 % yellow and green ambiguity cross...
 fill([piclim(1)-piclim(2) piclim(2)-piclim(1) piclim(2)-piclim(1) piclim(1)-piclim(2)],thick*[piclim(3)-piclim(4) piclim(3)-piclim(4) piclim(4)-piclim(3) piclim(4)-piclim(3)],'y');
 plot([piclim(1)-piclim(2) piclim(2)-piclim(1) piclim(2)-piclim(1) piclim(1)-piclim(2) piclim(1)-piclim(2)],thick*[piclim(3)-piclim(4) piclim(3)-piclim(4) piclim(4)-piclim(3) piclim(4)-piclim(3) piclim(3)-piclim(4)],[clr '-']);
 fill(thick*[piclim(1)-piclim(2) piclim(2)-piclim(1) piclim(2)-piclim(1) piclim(1)-piclim(2)],[piclim(3)-piclim(4) piclim(3)-piclim(4) piclim(4)-piclim(3) piclim(4)-piclim(3)],'y');
 plot(thick*[piclim(1)-piclim(2) piclim(2)-piclim(1) piclim(2)-piclim(1) piclim(1)-piclim(2) piclim(1)-piclim(2)],[piclim(3)-piclim(4) piclim(3)-piclim(4) piclim(4)-piclim(3) piclim(4)-piclim(3) piclim(3)-piclim(4)],[clr '-']);
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
  signcellred=signell(cell)*signell(1);
  signnoncellred=signell(noncell)*signell(1); 
  amosign=signnoncellred.*[prod(signcellred.^(spram(:,1)')),prod(signcellred.^(spram(:,2)'))];  
  % adjust the plot signs, according to whether 
  % you're plotting just the positive part, 
  % the unfolded real part, or the folded real part...
  plotsign=flip.*amosign; if fold==2 plotsign=[1 1]; end; 
 
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
    logabscellred1=logabsell1(cell)-logabsell1(1); 
     logabscellred2=logabsell2(cell)-logabsell2(1); 
    logabsnoncellred1=logabsell1(noncell)-logabsell1(1); 
     logabsnoncellred2=logabsell2(noncell)-logabsell2(1); 
 
    amoabs1=logabsnoncellred1+logabscellred1*pram; 
     amoabs2=logabsnoncellred2+logabscellred2*pram;  
    if ((amoabs1(1)>=piclim(1))&(amoabs1(1)<=piclim(2))&(amoabs1(2)>=piclim(3))&(amoabs1(2)<=piclim(4)))|((amoabs2(1)>=piclim(1))&(amoabs2(1)<=piclim(2))&(amoabs2(2)>=piclim(3))&(amoabs2(2)<=piclim(4)))  
     plot(plotsign(1)*([amoabs1(1) amoabs2(1)]-piclim(1)),plotsign(2)*([amoabs1(2) amoabs2(2)]-piclim(3)),[clr '-']); 
    end; % if you're in the plot window...
   end; % for j...
  end; % if the signs are correct...
  % adjust the window, according to whether
  % you're plotting just the positive part,
  % the unfolded real part, or the folded real part...
  if (fold==0)|(fold==2)
   % axis([0 piclim(2)-piclim(1) 0 piclim(4)-piclim(3)]); % change back!!!
  else % fold had better be 1 at this point...
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

axis equal; 

