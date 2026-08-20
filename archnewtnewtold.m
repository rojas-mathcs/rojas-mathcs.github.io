% univariate experiment for ARCHIMEDEAN amoeba thm! (april 24, 2004, 0:06)
%
% the numbers output spit out closely 
% predict the norms and phases! now the question is the deviation 
% from the predictions...
%
% july 22, 2004: the quality of the prediction appears to depend on 
%                (a) the ``conditioning'' of ArchNewt (are the 
%                    vertex angles nastily close to pi?) and 
%                (b) whether the maximal number of vertices 
%                    possible for the lower hull all appear
%
% NOTE: this program randomly picks coefficients, 
%       then sees how Newton iteration behaves, starting from 
%       the polyhedral predictions... 

% these are the exponents, written in descending order...
exps=[35,22,18,10,5,0];
% exps=[6 5 4 3 2 1 0];
% m = # of monomial terms = size of the above array... 
m=size(exps,2);
% d = degree (generically) of our random m-nomial = largest elemement of
%     the array exps
d=max(exps);
% these are the weights for the random gaussians used for the coefficients
% vars=[binom(exps(1),exps(:))]; 
vars=ones(1,m);
% how many random m-tuples of phases will we pick? 
% = # of iterations of newton...
num=50; 
% how many examples do you want?
exa=10;
% how well-conditioned do you demand the lower hull to be? 
% tol := # of edges * minimum difference between successive inner normals...
% (larger means better conditioned => closer match between guess and true) 
% tol=.3*pi;  
tol=.2*pi;  
% must insert blanks into array of coefficients so that
% matlab knows what to do with the missing monomial terms...
blanks=(exps(1:(m-1))-exps(2:m))-1;

tot=1;
% now generate random polynomials by randomly varying the 
% amplitudes of the coefficients...
while tot<=exa 

 % define the random coefficients...
 cuffs=(rand(1,m)+sqrt(-1)*rand(1,m))*diag([sqrt(vars(:))])';
 x=exps(:); y=-log(abs(cuffs))';
 % then get the Archimedean Newton polygon and its lower hull...
 k=convhull(x,y); mm=size(k,1);
 % isolate index of leftmost point
 [blah,indmin]=min(x(k(:)));
 % then rotate/reindex...
 lowerind=k(1+mod(indmin+[-1:(mm-3)],mm-1))';

% % don't bother plotting unless the lower hull is ``full''... 
% if (mm-1)==m 
 % why not allow any number of lower edges... 
 if 1==1

  % predict ord(roots):=log(roots)... 
  % (or is ord really defined as -log? check!) 
  ord=zeros(1,mm-2);
  for j=1:(mm-2)
   ord(j)=(y(lowerind(j+1))-y(lowerind(j)))/(x(lowerind(j+1))-x(lowerind(j)));
  end; 
 
  % ...and skip if lower hull is ``ill-conditioning''... 
  minangle=min([atan(ord([2:(mm-2)]))-atan(ord([1:(mm-3)]))]); 
  if minangle>tol/(mm-2)    

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
                                                                                
   % count of plotted polynomial roots goes up...
   tot=tot+1;
   % plot archimedean newton polygon and its lower hull...
   figure(2);
   plot(x(k(:)),y(k(:)),'k'); plot(x(lowerind(:)),y(lowerind(:)),'r');
   % for archnewts plot, size the axes, label them, and put on supertitle... 
   axis tight; axis equal;
   xlabel('Spectrum (or Support or Exponents)'); ylabel('-log(abs(coeff))');
   title(sprintf('%d-Edged Lower Hull of Archimedean Newton Polygon of a Random %d-nomial',mm-2,m));
  
   % define the lower inner normals in the form (ord,1),
   % then use ord as a guess for the absolute values of the roots...
   % then use the length of the edges, and adjacent coefficient ratios,
   % to make guesses (phi) for the phases of the roots...
   phi=zeros(2,mm-2); x0=[];
   for j=1:(mm-2)
    phi(2,j)=x(lowerind(j+1))-x(lowerind(j)); 
    phi(1,j)=angle(-cuffs(lowerind(j+1))/cuffs(lowerind(j)))/phi(2,j); 
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
 
   % plot the roots (according to matlab) of our m-nomial...
   plot(mod(pi+real(z),2*pi)-pi,imag(z),'b.'); 
   % now let's try to correctly plot the polyhedral predictions...
   for j=1:(mm-2)
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
   axis tight; title(sprintf('Roots of a Random %d-nomial with %d-Edged Lower Hull\n (BLUE = True Roots Differing From Polyhedral Prediction)\n (CYAN = %d Newton Steps From Polyhedral Prediction)',m,mm-2,num));
   % label the axes
   xlabel('Im(Log(Roots))'); ylabel('Re(Log(Roots))');
   if tot<=exa 
    pause;  
   end;
  end; % if minangle...
 end; % if mm-1==m...
end; % for i...

