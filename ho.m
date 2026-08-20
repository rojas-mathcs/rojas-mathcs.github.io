% univariate polynomial homotopy experiment: (june 21, 2004)
% copyright 2004, J. Maurice Rojas 
%
% This code fixes a support set and then illustrates 
% the well-known basic homotopy algorithm for solving 
% a univariate polynomial with given support. The homotopy  
% starts from the roots of x^d-1, where d is the degree. 
% 
% This program actually does this several times, awaiting 
% a key-stroke from the user to go on to the next example. 
% 

% these are the exponents, written in descending order...
exps=[22,10,8,5,0]; 
% m = # of monomial terms = size of the above array...
m=size(exps,2); 
% d = degree (generically) of our random m-nomial = largest elemement of
%     the array exps
d=max(exps);
% these are the weights for the random gaussians used for the coefficients
vars=[binom(exps(1),exps(:))]; 
% must insert blanks into array of coefficients so that
% matlab knows what to do with the missing monomial terms...
blanks=(exps(1:(m-1))-exps(2:m))-1;

% the coefficients of our simple start binomial...
doffs=[1,zeros(1,d-1),1];
% this is the number of examples we'll see...
num=10; 
% this is the number of steps used in each homotopy...
steps=100;

% the main loop: i indicates ith homotopy example...
for i=1:num

% erase figure 1 and get ready for plotting...
figure(1);
hold off;
subplot(1,1,1);
newplot;
hold on;

% erase figure 2 and get ready for plotting...
figure(2);
hold off;
subplot(1,1,1);
newplot;
hold on;

% initial guesses for the roots of our polynomial 
% are just the dth roots of unity...
% xi = the ith guess
% note that x0 and xi are ARRAYS consisting of all d guesses...
 x0=roots(doffs); xi=x0;  
% the coefficients of f --- the polynomial whose roots we want to find...
 cuffs=(randn(1,m)+sqrt(-1)*randn(1,m))*diag([sqrt(vars(:))]);  
% now successively insert zeroes so that matlab knows what to 
% do with the missing monomials...
 fcoeffs=[cuffs(1)]; 
 for j=1:(m-1)
  fcoeffs=[fcoeffs,zeros(1,blanks(j)),cuffs(j+1)];
 end; 

% now start the deformation of x^d - 1 into f...
 for j=1:(steps-1) 
  % the step size...
  t=j/steps; 
  % the coefficients of the jth deformation 
  % are just a convex linear combination of the coefficients of 
  % x^d-1 and f...
  coeffs=(1-t)*doffs+t*fcoeffs; 
  % evaluate our jth deformation polynomial 
  % via horner's rule, which is a bit more numerically stable 
  % than summing high powers...
  fxi=coeffs(1)*ones(d,1); 
  fpxi=d*coeffs(1)*ones(d,1); 
  for k=1:(d-1)
   fxi=(xi.*fxi)+coeffs(k+1)*ones(d,1); 
   fpxi=(xi.*fpxi)+(d-k)*coeffs(k+1)*ones(d,1);
  end;
  fxi=(xi.*fxi)+coeffs(d+1)*ones(d,1);  
  % update our root guesses by one iteration of newton, 
  % and plot the jth GUESSES in cyan in figure 1...
  xi=xi-(fxi./fpxi); figure(1); plot(xi,'c.');  
  % get the TRUE ROOTS of the jth deformation and plot them in yellow 
  % in figure 2...
  r=roots(coeffs); figure(2); plot(r,'y.');
 end;
 figure (1);
 % plot the roots of f in blue...
 r=roots(fcoeffs); plot(r,'b*'); 
 % plot the roots of x^d-1 in red...
 % make sure axes aren't stretched...
 r=roots(doffs); plot(r,'r*');  axis equal; 
 % make supertitle...
 title(sprintf('%d Step Homotopy Between Degree %d Random %d-nomial and x^{%d}-1',steps,d,m,d));
 % label axes...
 xlabel('Re(Root)');
 ylabel('Im(Root)'); 
 figure (2); 
 % plot the roots of f in blue...
 r=roots(fcoeffs); plot(r,'b*'); 
 % plot the roots of x^d-1 in red...
 % make sure axes aren't stretched...
 r=roots(doffs); plot(r,'r*');  axis equal; 
 % make supertitle...
 title(sprintf('%d Step Homotopy Between Degree %d Random %d-nomial and x^{%d}-1',steps,d,m,d));
 % label axes...
 xlabel('Re(Root)');
 ylabel('Im(Root)'); 
 % if we haven't already done num many examples, 
 % wait for a keystroke to do another...
 if num>1 pause; end;  
end; 
