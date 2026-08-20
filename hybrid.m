% sets output to high precision
format('long');
% clears variable assignments and command window
clear;
clc;
% Set degree d of polynomial
d = 20;
% Set coefficient row vector C of polynomial
C = randn(1,d+1);
for j=1:d+1
    if rand(1)>0.75
        C(j)=C(j);%.*sqrt(mfun('binomial',d,j-1));
    else
        C(j) = 0;
    end
end
display(C);

% HSS algorithm

tic;
% s is the humber of circles needed
s = ceil(.26632*log(d));
% N is the number of points on each circle
N = ceil(8.32547*d*log(d));
% calculate the starting points for Newton iteration, store them in R
R = [];
for v=1:s
    for j=0:N-1
        R = [R; exp((2*pi*i*j)./N)*(1+sqrt(2)).*(((d-1)./d).^((2.*v-1)./(4.*s)))];
    end
end
% Newton iteration using the roots in the column vector R as start points
% nmax is the maximum number of iterations
nmax=400;
% j counts the numbers of iterations
j=0;
%iteration loop
while j<nmax
    R = [R, R(:,j+1)-(polyval(C,R(:,j+1))./polyval(polyder(C),R(:,j+1)))];
    j = j+1;
end

% eliminates redundant roots from array of calculated roots
resultList1 = (10^10).*R(:,j+1);
resultList1 = round(resultList1);
resultList1 = resultList1./(10^10);
resultList1 = unique(resultList1);
display(resultList1);
time = toc

% Rojas algorithm

tic
% eliminate zero-valued coefficients; calculate Archimedean Newton polygon
% M contains the first coordinates, and N contains the second coordinates
M = [];
N = [];
for j = 1:d+1
    if C(j) ~= 0
        M = [M,j-1];
        N = [N, -log(abs(C(j)))];
    end
end
        
% computes column vector E of convex hull vertex indices in M
E=convhull(M,N);
[rows,dummy] = size(E);
[dummy, rightmost] = size(M);

% crop last element of E, since it is a repeat of the first element
E=E(1:rows-1);

% reduces column vector E to LOWER convex hull vertex indices in M
for a=1:rows-1
    if E(a)==1
        for j=a:rows-1
            if E(j)==rightmost
                E=E(a:j);
                k=j-a+1;
                break;
            end
        end
        break;
    elseif E(a)==rightmost
        for j=a:rows-1
            if E(j)==1
                E=[E(j:rows-1);E(1:a)];
                k=rows-j+a;
                break;
            end
        end
        break;
    end
end

% solve for the complex roots, place into column vector R
R=[];
for a=1:k-1
    R = [R;roots([1,zeros(1,M(E(a+1))-M(E(a))-1),C(M(E(a+1))+1)/C(M(E(a))+1)])];
end

% Newton iteration using the roots in the column vector R as start points
% nmax is the maximum number of iterations
nmax=400;
% j counts the numbers of iterations
j=0;
%iteration loop
while j<nmax
    R = [R, R(:,j+1)-(polyval(C,R(:,j+1))./polyval(polyder(C),R(:,j+1)))];
    j = j+1;
end
resultList2 = (10^10).*R(:,j+1);
resultList2 = round(resultList2);
resultList2 = resultList2./(10^10);
resultList2 = unique(resultList2);
display(resultList2);
time = toc

% resets output to low precision
format('short');
