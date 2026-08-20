% univariate polynomial root music experiment!: (june 21, 2004)
% copyright 2004, J. Maurice Rojas...
% 
% This program computes the roots of many random polynomials, with 
% fixed monomial term structure, and gives 3 interesting plots: 
% figure 1: the roots themselves (a pretty wheel-like picture emerges) 
% figure 2: the log of the roots (the wheel is unwound and the spokes are 
%                                 clarified) 
% figure 3: a histogram based on figure 2, showing the density of the random 
%           roots.
% this illustrates the ``harmonious'' fundamental theorem of algebra: 
% the fact that norms and phases of the roots are determined by 
% a collection of binomials determined by the lower hull of the 
% newton polygon. in particular, one sees that the norms of the roots occur in 
% m regions (where m is the number of monoimal terms) and the phases 
% of the roots depend on the pairs of coefficients defined by a lower 
% edges of the Archimedean Newton polygon.  
%
% As an added bonus, the histogram in figure 3 is interpreted as 
% the frequency spectrum of an audio signal, and the resulting 
% audio data is placed in two arrays: sigadd and sigmult. 
% These can be played with the command sound(sigadd) or sound(sigmult). 
% 
% The underlying probability measure weighs the middle coefficients 
% more greatly (more on this later in class). 

% get figure 1 ready...
figure(1); 
hold off; 
subplot(1,1,1);
newplot;
hold on; 

% get figure 2 ready...
figure(2); 
hold off; 
subplot(1,1,1);
newplot;
hold on; 

% get figure 3 ready...
figure(3); 
hold off; 
subplot(1,1,1);
newplot;
hold on; 

% these are the exponents, written in descending order...
exps=[31,23,22,20,0];
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
% number of random m-nomials to try
num=200; 
% number of samples for the audio signal in the time dimension...
tym=100;
% how many different frequencies are we using for our audio signal?
freq=72; 
numnotes=72; 
% lowest frequency...
flow=440/(2^(34/12)); % 34 half steps below concert A
% highest frequency...
fhi=440*(2^(37/12)); % 37 half steps above concert A 

% get the audio signal array ready...
spec=zeros(tym,freq);
% get the array of m-nomial roots ready...
z=zeros(exps(1),num);

% make sure the next plot goes to figure 1...
figure(1); 
% the main outer loop: do something with num many 
% random m-nomials...
for i=1:num
 % generate the coefficients we need, as POSITIVE weighed gaussians...
 cuffs=abs(randn(1,m)*diag([sqrt(vars(:))])); 
 % insert the missing 0 coefficients successively...
 coeffs=[cuffs(1)];
 for j=1:(m-1)
  coeffs=[coeffs,zeros(1,blanks(j)),cuffs(j+1)];
 end;
 % solve for the roots of our random m-nomial 
 r=roots(coeffs); 
 % plot the roots in figure 1 
 figure(1); plot(r,'.'); 
 % the ith column of z is the vector consisting of the lots of the roots 
 z(:,i)=sqrt(-1)*log(r);
 % plot log(roots) in figure 2 
 figure(2); plot(z(:,i),'.');
end; 
% make sure the axes of figure 1 have equal proportions
figure(1); axis equal;  
% put the supertitle on the plot 
title(sprintf('Roots of %d Random %-nomials',num,m));
% label the axes 
xlabel('Re(Root)');
ylabel('Im(Root)'); 
% make sure the axes of figure 2 have equal proportions
figure(2); axis equal;  
% put the supertitle on the plot 
title(sprintf('Sqrt(-1)*Log(Roots) of %d Random %d-nomials',num,m));
% label the axes 
xlabel('Re(Root)');
ylabel('Im(Root)'); 

% find where the audio signal should start/end by 
% finding the largest and smallest root norms...
mintime=min(min(imag(z)));
maxtime=max(max(imag(z))); 

% this creates the histogram: if an element of z falls 
% into a certain box, then a certain element of the 
% histogram/audio signal array (spec) is increased by 1. 
% After running through the outer loop, this counts the number of roots 
% falling in all the boxes...
for i=1:num
 for j=1:exps(1)
  coords=1+floor([real(z(j,i))+pi,imag(z(j,i))-mintime]*[(freq-1)/(2*pi) 0 ; 0 (tym-1)/(maxtime-mintime)]);  
  spec(coords(2),coords(1))=spec(coords(2),coords(1))+1;
 end; 
end;

% find the sums of the elemnts in the rows of spec...
mxs=spec*ones(freq,1); 
% find the largest row sum 
mx=max(mxs); 

figure(3)
% plot the histogram as an interpolated colored surface...
h=surfl(spec); 
shading interp;
% make sure the axes are scaled so that as much as possible of the 
% surface is shown...
axis tight;   
% insert supertitle
title(sprintf('Rescaled Histogram of Sqrt(-1)*Log(Roots) for %d Random %d-nomials',num,m));
% label the axes...
xlabel('Arg(Root) renormalized');
ylabel('Log(Abs(Root)) renormalized');
% NOTE: don't forget to click on the little curled arrow 
%       on the top right of figure 3 so you can click 
%       and drag to make rotations of the 3D figure! 

% the rows of spec will now be turned into audio blips...
basenote=flow*[0:1023]*pi/4096; % this gives 1/8 second per row...
% this is the array of allowed frequencies, piano-style...
notesmult=exp([0:log(fhi/flow)/(freq-1):log(fhi/flow)])'; 
% this is the array of allowed frequencies, signal processing style...
notesadd=[1:(fhi/flow)/(freq+1):fhi/flow]';
% compute the necessary sine waves constituting our audio signals 
harmadd=sin(notesadd*basenote); 
harmmult=sin(notesmult*basenote);

% normalize our audio signal 
spec=spec/mx;
% compute the necessary linear combinations of sine waves to 
% get our audio signals! 
stackedsigadd=(spec*harmadd)'; 
stackedsigmult=(spec*harmmult)'; 
sigadd=stackedsigadd(:);
sigmult=stackedsigmult(:);

