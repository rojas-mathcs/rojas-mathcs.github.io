% Copyright 2008, J. Maurice Rojas
% modified archnewtnewt to give approximations of NONzero roots 
% of a sparse polynomial...
% 
% usage: [xi,ord,phi]=archnewtnewt(exps,num,cuffs); 
% where... 
% exps   = array of exponents, in decreasing order...
% num    = # of times to do newton...
% cuffs  = coefficients for sparse polynomial in DESCENDING exponent order...  
%
% exps=[35,31,18,14,8,0]; num=50; cuffs=[1 5 10 10 5 1]; 
% should give a non-trivial output...
% note: the output is in the following format: 
%         xi = array of roots, after applying num iterations of 
%               newton to the lower binomial guesses...
%        ord = [log|norm of cluster #1|  log|norm of cluster #2| ... ] 
%        phi = [ argumentcluster#1 argumentcluster#2 ... ] 
function [xi,ord,phi]=archnewtnewt(exps,num,cuffs) 

% eliminate repeats in exponents...
mm=size(exps,2); y=unique(exps); m=size(y,2); exps=y(m:-1:1); 
% m = # of monomial terms = size of the above array... 

if mm>m 
 disp(sprintf('Had to eliminate some repeats in your exponents...\n')); 
end; 

% d = degree (generically) of our random m-nomial = largest elemement of
%     the array exps
d=exps(1);

% let's try to detect colinearity of 3 consecutive lifted support points... 
tol=1e-7; 

% must insert blanks into array of coefficients so that
% matlab knows what to do with the missing monomial terms...
blanks=(exps(1:(m-1))-exps(2:m))-1; 
if exps(m)>0 % fix case of zero constant term...
 d=d-exps(m); 
 exps=exps-exps(m); 
end;  

 % define the random coefficients...
 x=exps(:); y=-log(abs(cuffs))';

 if m>2 
  for i=1:(m-2) 
   if abs(det([1 1 1; x(i) x(i+1) x(i+2); y(i) y(i+1) y(i+2)]))<tol 
    xi=[]; ord=[]; phi=[]; 
    disp(sprintf('Three lifted support points are colinear! Not good!\n')); 
    return; 
   end; 
  end;  
 end; 
 
 % then get the Archimedean Newton polygon and its lower hull...
 k=convhull(x,y); mm=size(k,1); 
 % isolate index of leftmost point
 [blah,indmin]=min(x(k(:)));
 % then rotate/reindex...
 lowerind=k(1+mod(indmin+[-1:(mm-3)],mm-1))'; 
 [blah,mm]=min(lowerind); % find position of 1! i.e., mm = # of lower verts! 

  % predict ord(roots):=log(roots)... 
  % (or is ord really defined as -log? check!) 
  ord=zeros(1,mm-1);
  for j=1:(mm-1)
   ord(j)=(y(lowerind(j+1))-y(lowerind(j)))/(x(lowerind(j+1))-x(lowerind(j)));
  end; 
 
   % delete the following 11 lines if you want to remove plotting...
   figure(1);
   hold off;
   subplot(1,1,1);
   newplot;
   hold on;

   figure(2);
   hold off;
   subplot(1,1,1);
   newplot;
   hold on;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % if you want to remove plotting, also delete all lines 
   % that say ``plot'', ``figure'', ``axis'', ``title'', or ``xlabel''...

   % plot archimedean newton polygon and its lower hull...
   figure(2);
   plot(x(k(:)),y(k(:)),'k'); plot(x(lowerind(1:mm)),y(lowerind(1:mm)),'r-');
   % for archnewts plot, size the axes, label them, and put on supertitle... 
   axis tight; axis equal;
   xlabel('Spectrum (or Support or Exponents)'); ylabel('-log(abs(coeff))');
   title(sprintf('%d-Edged Lower Hull of Archimedean Newton Polygon of a Random %d-nomial',mm-1,m));
  
   % define the lower inner normals in the form (ord,1),
   % then use ord as a guess for the absolute values of the roots...
   % then use the length of the edges, and adjacent coefficient ratios,
   % to make guesses (phi) for the phases of the roots...
   phi=zeros(2,mm-1); x0=[];
   for j=1:(mm-1)
    phi(2,j)=x(lowerind(j+1))-x(lowerind(j));  % horizontal length of 
                                               % lower edge!!! 
    phi(1,j)=angle(-cuffs(lowerind(j+1))/cuffs(lowerind(j)))/phi(2,j); 
           % vertical length of lower edge! 
    for l=0:(phi(2,j)-1)
     x0=[x0;exp(-sqrt(-1)*(phi(1,j)+l*2*pi/phi(2,j))+ord(j))]; 
    end;
   end; 

   % insert the missing 0 coefficients successively...
   coeffs=[cuffs(1)];
   for j=1:(m-1)
    coeffs=[coeffs,zeros(1,blanks(j)),cuffs(j+1)];
   end;
   z=sqrt(-1)*log(roots(coeffs));

   % now try to do a little newton iteration, employing 
   % horner's rule along the way...
   xi=x0; figure(1); 
   for j=1:num 
    fxi=coeffs(1)*ones(d,1);
    fpxi=d*coeffs(1)*ones(d,1);
 
    for k=1:(d-1)
     fxi=(xi.*fxi)+coeffs(k+1)*ones(d,1);
     fpxi=(xi.*fpxi)+(d-k)*coeffs(k+1)*ones(d,1);
    end;
    fxi=(xi.*fxi)+coeffs(d+1)*ones(d,1);
    % update our root guesses by one iteration of newton,
    % and plot the jth GUESSES in cyan in figure 1...
    xi=xi-(fxi./fpxi); plot(sqrt(-1)*log(xi),'c.');
   end; 

   % plot ostrowski's strips!!!
   nu=0; % ostrowski's parameter starts at right endpoint of 1st edge!!! 
   for j=1:(mm-1)
    nu=nu+phi(2,j);  
    plot([-pi pi pi -pi],ord(j)+[log(1-.5^(1/nu))+[0 0] [0 0]-log(1-.5^(1/(d-nu+1)))],'m-'); 
    fill([-pi pi pi -pi],ord(j)+[log(1-.5^(1/nu))+[0 0] [0 0]-log(1-.5^(1/(d-nu+1)))],'m'); 
   end; 
 
   % plot the roots (according to matlab) of our m-nomial...
   plot(mod(pi+real(z),2*pi)-pi,imag(z),'b.'); 
   % now let's try to correctly plot the polyhedral predictions...
   for j=1:(mm-1)
    if phi(1,j)>=0
     for l=0:(phi(2,j)-1)
      plot(mod(pi+phi(1,j)+l*2*pi/phi(2,j),2*pi)-pi,ord(j),'r.');
     end;
    else
     for l=0:(phi(2,j)-1)
      plot(mod(3*pi+phi(1,j)+l*2*pi/phi(2,j),2*pi)-pi,ord(j),'r.');
     end;
    end;
   end; 
    
   % final prettying up of the roots plot... 
   axis tight; title(sprintf('Roots of a Random %d-nomial with %d-Edged Lower Hull\n BLUE = True Roots , RED= Polyhedral Prediction\n CYAN = %d Newton Steps From Polyhedral Prediction',m,mm-1,num));
   % label the axes
   xlabel('Im(Log(Roots))'); ylabel('Re(Log(Roots))');

