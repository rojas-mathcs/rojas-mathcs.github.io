% univariate experiment for ARCHIMEDEAN amoeba thm! (april 24, 2004, 0:06)
%
% the numbers output closely 
% predict the phases! now the question is the deviation 
% from the predictions...
%
% july 22, 2004: the quality of the prediction appears to depend on 
%                (a) the ``conditioning'' of ArchNewt (are the 
%                    vertex angles nastily close to pi?) and 
%                (b) whether the maximal number of vertices 
%                    possible for the lower hull all appear
%
% NOTE: this program fixes one of coefficient phases, 
%       THEN varies the amplitudes randomly...
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

% these are the exponents, written in descending order...
exps=[35,31,18,14,8,0];
% m = # of monomial terms = size of the above array... 
m=size(exps,2);
% d = degree (generically) of our random m-nomial = largest elemement of
%     the array exps
d=max(exps);
% these are the weights for the random gaussians used for the coefficients
vars=[binom(exps(1),exps(:))];
% vars=ones(1,m);
% how many random m-tuples of phases will we pick? 
% = number of random m-nomials to try
num=50; 
% must insert blanks into array of coefficients so that
% matlab knows what to do with the missing monomial terms...
blanks=(exps(1:(m-1))-exps(2:m))-1;

% generate the coefficient absolutes we'll fix (positive Gaussians)...
phases=exp(2*pi*sqrt(-1)*randn(1,m))'; 


% now generate random polynomials by randomly varying the 
% amplitudes of the coefficients...
for i=1:num

 % define archimedean lifting of the spectrum...
 amps=abs(randn(1,m)*diag([sqrt(vars(:))]))';
 x=exps(:); y=-log(amps);
 % then get the Archimedean Newton polygon and its lower hull...
 k=convhull(x,y); mm=size(k,1);
 % isolate index of leftmost point
 [blah,indmin]=min(x(k(:)));
 % then rotate/reindex...
 lowerind=k(1+mod(indmin+[-1:(mm-3)],mm-1))';

 % plot archimedean newton polygon and its lower hull...
 figure(2);
 plot(x(k(:)),y(k(:)),'k'); plot(x(lowerind(:)),y(lowerind(:)),'r');

 duffs=amps.*phases;  
 % insert the missing 0 coefficients successively...
 coeffs=[duffs(1)];
 for j=1:(m-1)
  coeffs=[coeffs,zeros(1,blanks(j)),duffs(j+1)];
 end;
 z=roots(coeffs); figure(1); plot(sqrt(-1)*log(z),'b.'); 

 % define the lower inner normals in the form (ord,1),
 % then use ord as a guess for the absolute values of the roots...
 % then use the length of the edges, and adjacent coefficient ratios,
 % to make guesses (phi) for the phases of the roots...
 ord=zeros(1,mm-2);
 phi=zeros(2,mm-2);
 for j=1:(mm-2)
  ord(j)=(y(lowerind(j+1))-y(lowerind(j)))/(x(lowerind(j+1))-x(lowerind(j)));
  phi(2,j)=x(lowerind(j+1))-x(lowerind(j)); 
  phi(1,j)=angle(-phases(lowerind(j+1))/phases(lowerind(j)))/phi(2,j); 
  for l=0:(phi(2,j)-1)
   plot(mod(pi+phi(1,j)+l*2*pi/phi(2,j),2*pi)-pi,ord(j),'r.');
  end;
 end;
end; 

% put the supertitle on the root plot
figure(1); 
axis tight; title(sprintf('Roots of %d Random %d-nomials with fixed coefficient phases\n (BLUE = True Roots Differing from Polyhedral Prediction)',num,m));
% label the axes
xlabel('Im(Log(Roots))'); ylabel('Re(Log(Roots))');

% for archnewts plot, size the axes, label them, and put on supertitle... 
figure(2); axis tight; axis equal;
xlabel('Spectrum (or Support or Exponents)'); ylabel('-log(abs(coeff))');
title(sprintf('Lower Hulls of Archimedean Newton Polygons of %d Random %d-nomials',num,m));
