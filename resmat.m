% Program Description: This program helps to find the mixed subdivision
%                      of three polynomials in two unknowns. It also 
%                      form a toric resultant matrix corresponding
%                      to a polynomial system.
%
% This program was written by Christopher Wong Siu Cheung 
% for a final year project directed by J. Maurice Rojas. 
% Last modified on 2000/04/28.
%
% Points to note:
% 1. If you need to solve 2 by 2 polynomial system,
%    an input format you can use:
%
%    f1 = [deg(1st term in f1) deg(2nd term in f1)...]
%    f2 = [deg(1st term in f2) deg(2nd term in d2)...]
%    f3 = [0 0 1; 0 1 0]
%
% 2. Output File "res.m" : The toric resultant matrix of corresponding 
%                          to the input polynomial system.
%

function resmat32(P, Q, R)
global mcptlst;
global umcptlst;
global eps;
global cellptlst;
global mixed;
global lattice;
global lattice_info;

TRUE = 1;
FALSE = 0;

color{1} = [0 0 1];
color{2} = [1 0 0];
color{3} = [0 1 0];
color{4} = [1 1 0.6];
color{5} = [0.8 1 1];
color{6} = [1 0.8 1];

% The allarea array stored the total area of unmixed and mixed cell
% (ie. allarea = [blue red green yellow cyan pink])
allarea = [0 0 0 0 0 0];

% Set up an array for mixed cell's lattice point list
mcptlst = [];

% Set up a data structure for storing the cell point list
% (ie. cellptlst = {{blue}, {red}, {green}, {yellow}, {cyan}, {pink}}
cellptlst = {[],[],[],[],[],[]}

% Set up a data structure for storing the lattice point
% (ie. lattice = [x0 x1 x2 ... 
%                 y0 y1 y2 ... ] )
%                  ^lattice point in cell 0
%
% lattice point = integer points inside the cell
lattice = []

% Set up a data structure for storing the information of lattice point
% (ie. lattice_info = [x0 x1 x2 ...
%                      y0 y1 y2 ...
%                      c0 c1 c2 ... ])
% (x,y) is coordinate of cell c
lattice_info = []

[hullP, nvpolyP, edgesP, polyP] = analysis(P);
[hullQ, nvpolyQ, edgesQ, polyQ] = analysis(Q);
[hullR, nvpolyR, edgesR, polyR] = analysis(R);
poly = {[polyP], [polyQ], [polyR]};
nvpoly = [nvpolyP nvpolyQ nvpolyR]

% FOR UNMIXED AREA
close ALL;
hold on;
[cutp1, cutp2, cutq1, cutq2, cutr1, cutr2] = selectcut(polyP, polyQ, polyR);
cut = {[cutp1 cutp2], [cutq1 cutq2], [cutr1 cutr2]}
cutedge = [poly{1}(:,cut{1}(2))-poly{1}(:,cut{1}(1)),...
      poly{2}(:,cut{2}(2))-poly{2}(:,cut{2}(1)),...
      poly{3}(:,cut{3}(2))-poly{3}(:,cut{3}(1))]

if (det([cutedge(:,1), cutedge(:,2)])>0)
   seq = [1 2];
else
   seq = [2 1];
end

% Generally, 4 possible cases.
if (det([cutedge(:,seq(1)), cutedge(:,3)])>0)
   if (det([cutedge(:,3), cutedge(:,seq(2))])>0)
      % Case 1
      temp = cut{3}(1);
      cut{3}(1) = cut{3}(2);
      cut{3}(2) = temp;
      cutedge(:,3) = -cutedge(:,3);
      seq = [seq(1) seq(2) 3];
   else
      % Case 2
      temp = cut{seq(2)}(1);
      cut{seq(2)}(1) = cut{seq(2)}(2);
      cut{seq(2)}(2) = temp;
      cutedge(:,seq(2)) = -cutedge(:,seq(2));
      seq = [seq(1) 3 seq(2)];
   end
else
   if (det([cutedge(:,3), cutedge(:,seq(2))])>0)
      % Case 4
      temp = cut{seq(1)}(1);
      cut{seq(1)}(1) = cut{seq(1)}(2);
      cut{seq(1)}(2) = temp;
      cutedge(:,seq(1)) = -cutedge(:,seq(1));
      seq = [seq(1) 3 seq(2)];
   else
      % Case 3
      seq = [seq(1) seq(2) 3];
   end
end

% Calculate eps = [epsx epsy], where epsx and epsy are two real numbers between 0 and 1
% All the polygon will be shifted by x->x+epsx and y->y+eps
eps = rand(1,2);

% Testing data
% eps = [0.8214 0.4447];
% eps = [0 0];

% FORM ALL UNMIXED CELL, NOTE THAT TOTALLY 6 UNMIXED CELLS
% Note that the data structure of storing mixed and unmixed cell is different!
unmixed = {
   {[genindex(cut{seq(1)}(1), cut{seq(1)}(2), nvpoly(seq(1)))], [cut{seq(2)}(1)], [cut{seq(3)}(2)]},...
   {[genindex(cut{seq(1)}(2), cut{seq(1)}(1), nvpoly(seq(1)))], [cut{seq(2)}(2)], [cut{seq(3)}(1)]},...
   {[cut{seq(1)}(2)], [genindex(cut{seq(2)}(1), cut{seq(2)}(2), nvpoly(seq(2)))], [cut{seq(3)}(1)]},...
   {[cut{seq(1)}(1)], [genindex(cut{seq(2)}(2), cut{seq(2)}(1), nvpoly(seq(2)))], [cut{seq(3)}(2)]},...
   {[cut{seq(1)}(2)], [cut{seq(2)}(1)], [genindex(cut{seq(3)}(2), cut{seq(3)}(1), nvpoly(seq(3)))]},...
   {[cut{seq(1)}(1)], [cut{seq(2)}(2)], [genindex(cut{seq(3)}(1), cut{seq(3)}(2), nvpoly(seq(3)))]}
}

% Each cell data structure {[index in seq(1)], [index in seq(2)], [index in seq(3)]}
mixed = {};

% Fill the unmixed area for the all the unmixed cell list
unmixed_loop = {
   {[genindex(cut{seq(1)}(1), cut{seq(1)}(2), nvpoly(seq(1))) cut{seq(1)}(1)], [cut{seq(2)}(1)], [cut{seq(3)}(2)]},...
   {[genindex(cut{seq(1)}(2), cut{seq(1)}(1), nvpoly(seq(1))) cut{seq(1)}(2)], [cut{seq(2)}(2)], [cut{seq(3)}(1)]},...
   {[cut{seq(1)}(2)], [genindex(cut{seq(2)}(1), cut{seq(2)}(2), nvpoly(seq(2))) cut{seq(2)}(1)], [cut{seq(3)}(1)]},...
   {[cut{seq(1)}(1)], [genindex(cut{seq(2)}(2), cut{seq(2)}(1), nvpoly(seq(2))) cut{seq(2)}(2)], [cut{seq(3)}(2)]},...
   {[cut{seq(1)}(2)], [cut{seq(2)}(1)], [genindex(cut{seq(3)}(2), cut{seq(3)}(1), nvpoly(seq(3))) cut{seq(3)}(2)]},...
   {[cut{seq(1)}(1)], [cut{seq(2)}(2)], [genindex(cut{seq(3)}(1), cut{seq(3)}(2), nvpoly(seq(3))) cut{seq(3)}(1)]}
}
% Set up an array for unmixed cell's lattice point list
umcptlst = [];
for i=1:6
   xord = poly{seq(1)}(1,unmixed_loop{i}{1})+poly{seq(2)}(1,unmixed_loop{i}{2})+poly{seq(3)}(1,unmixed_loop{i}{3});
   yord = poly{seq(1)}(2,unmixed_loop{i}{1})+poly{seq(2)}(2,unmixed_loop{i}{2})+poly{seq(3)}(2,unmixed_loop{i}{3});
   fill(xord+eps(1),yord+eps(2),color{seq(ceil(i/2))});
   plot(xord+eps(1),yord+eps(2),ccolor(seq(ceil(i/2))));
   allarea(seq(ceil(i/2))) = allarea(seq(ceil(i/2))) + polyarea(xord, yord);
   
   xord = poly{seq(1)}(1,unmixed{i}{1})+poly{seq(2)}(1,unmixed{i}{2})+poly{seq(3)}(1,unmixed{i}{3});
   yord = poly{seq(1)}(2,unmixed{i}{1})+poly{seq(2)}(2,unmixed{i}{2})+poly{seq(3)}(2,unmixed{i}{3});
   intpt = polyintpt([xord+eps(1);yord+eps(2)],0,seq(ceil(i/2)));
   
   if (intpt>0)
      formcell = zeros(1,3)
      formcell(seq(1:3)) = [size(unmixed{i}{1},2),size(unmixed{i}{2},2),size(unmixed{i}{3},2)]
      [u,v] = min(formcell)
      xord = poly{v}(1,unmixed{i}{find(seq==v)})
      yord = poly{v}(2,unmixed{i}{find(seq==v)})
      li = [xord, yord, v]'*ones(1,intpt)
      lattice_info = [lattice_info,li]
   end
end

% COMPUTE MIXED AREA
edges = {[edgesP], [edgesQ], [edgesR]}

% We already have 3 mixed cells initially in the middle formed by the 6 unmixed cells
% Fill those 3 mixed cells
tempmc{seq(1)} = [cut{seq(1)}(1)]
tempmc{seq(2)} = [cut{seq(2)}(1:2)]
tempmc{seq(3)} = [cut{seq(3)}(1:2)]
allarea = colorcell(tempmc, poly, seq, color, seq(1), seq(2), seq(3), allarea, 1);
mixed = {mixed{:},{tempmc{:}}}

tempmc{seq(1)} = [cut{seq(1)}(1:2)]
tempmc{seq(2)} = [cut{seq(2)}(1)]
tempmc{seq(3)} = [cut{seq(3)}(1:2)]
allarea = colorcell(tempmc, poly, seq, color, seq(2), seq(3), seq(1), allarea, 1);
mixed = {mixed{:},{tempmc{:}}}

tempmc{seq(1)} = [cut{seq(1)}(1:2)]
tempmc{seq(2)} = [cut{seq(2)}(1:2)]
tempmc{seq(3)} = [cut{seq(3)}(1)]
allarea = colorcell(tempmc, poly, seq, color, seq(3), seq(1), seq(2), allarea, 1);
mixed = {mixed{:},{tempmc{:}}}

% Note that the first round is trial
trial = TRUE;
% Trial Loop
bcellseq = 3;
fcellseq = mod(bcellseq-1-1,3)+1;
rcellseq = mod(bcellseq-2-1,3)+1;
fstart = cut{seq(fcellseq)}(2);
bstart = cut{seq(bcellseq)}(2);
fstop = cut{seq(fcellseq)}(1);
bstop = cut{seq(bcellseq)}(1);
cutr = 1;
rrecord = cut{seq(rcellseq)}(1);
% Clear the variable "fstoparray"
clear fstoparray;
[fstop, bstop, brecord, fstoparray, allarea] = genmixedcell(trial,fcellseq,bcellseq,rcellseq,fstart,fstop,bstart,bstop,cutr,rrecord,color,poly,nvpoly,cut,cutedge,seq,edges,allarea);
% Flag 1
% Record bstop for the last mixed-cell use
keepbstop = bstop
% Swap the front cell and the remain cell and do mixed area computation again
temp = fcellseq;
fcellseq = rcellseq;
rcellseq = temp;
fstart = cut{seq(fcellseq)}(1);
rrecord = fstop;
rrecordarray = fstoparray;
allarea = special_genmixedcell(fcellseq,bcellseq,rcellseq,fstart,bstart,rrecord,rrecordarray,color,poly,nvpoly,cut,cutedge,seq,edges,allarea);

trial = FALSE;
% Flag 2
bcellseq = 2;
fcellseq = mod(bcellseq-1-1,3)+1;
rcellseq = mod(bcellseq-2-1,3)+1;
fstart = cut{seq(fcellseq)}(1);
bstart = cut{seq(bcellseq)}(1);
bstop = fstop;
fstop = cut{seq(fcellseq)}(2);
cutr = 2;
rrecord = brecord;
% Clear the variable "fstoparray"
clear fstoparray;
[fstop, bstop, brecord, fstoparray, allarea] = genmixedcell(trial,fcellseq,bcellseq,rcellseq,fstart,fstop,bstart,bstop,cutr,rrecord,color,poly,nvpoly,cut,cutedge,seq,edges, allarea);
% Flag 3
% Swap the front cell and the remain cell and do mixed area computation again
temp = fcellseq;
fcellseq = rcellseq;
rcellseq = temp;
fstart = cut{seq(fcellseq)}(2);
rrecord = fstop;
rrecordarray = fstoparray;
allarea = special_genmixedcell(fcellseq,bcellseq,rcellseq,fstart,bstart,rrecord,rrecordarray,color,poly,nvpoly,cut,cutedge,seq,edges,allarea);

% Flag 4
bcellseq = 1;
fcellseq = mod(bcellseq-1-1,3)+1;
rcellseq = mod(bcellseq-2-1,3)+1;
fstart = cut{seq(fcellseq)}(2);
bstart = cut{seq(bcellseq)}(2);
bstop = fstop;
fstop = cut{seq(fcellseq)}(1);
cutr = 1;
rrecord = brecord;
% Clear the variable "fstoparray"
clear fstoparray;
[fstop, bstop, brecord, fstoparray, allarea] = genmixedcell(trial,fcellseq,bcellseq,rcellseq,fstart,fstop,bstart,bstop,cutr,rrecord,color,poly,nvpoly,cut,cutedge,seq,edges, allarea);
% Flag 5
% Swap the front cell and the remain cell and do mixed area computation again
temp = fcellseq;
fcellseq = rcellseq;
rcellseq = temp;
fstart = cut{seq(fcellseq)}(1);
rrecord = fstop;
rrecordarray = fstoparray;
allarea = special_genmixedcell(fcellseq,bcellseq,rcellseq,fstart,bstart,rrecord,rrecordarray,color,poly,nvpoly,cut,cutedge,seq,edges,allarea);

% Flag 6
bcellseq = 3;
fcellseq = mod(bcellseq-1-1,3)+1;
rcellseq = mod(bcellseq-2-1,3)+1;
fstart = cut{seq(fcellseq)}(1);
bstart = cut{seq(bcellseq)}(1);
bstop = fstop;
fstop = cut{seq(fcellseq)}(2);
cutr = 2;
rrecord = brecord;
% Clear the variable "fstoparray"
clear fstoparray;
[fstop, bstop, brecord, fstoparray, allarea] = genmixedcell(trial,fcellseq,bcellseq,rcellseq,fstart,fstop,bstart,bstop,cutr,rrecord,color,poly,nvpoly,cut,cutedge,seq,edges,allarea);
% Flag 7
% Swap the front cell and the remain cell and do mixed area computation again
temp = fcellseq;
fcellseq = rcellseq;
rcellseq = temp;
fstart = cut{seq(fcellseq)}(2);
rrecord = fstop;
rrecordarray = fstoparray;
allarea = special_genmixedcell(fcellseq,bcellseq,rcellseq,fstart,bstart,rrecord,rrecordarray,color,poly,nvpoly,cut,cutedge,seq,edges,allarea);

% Flag 8
bcellseq = 2;
fcellseq = mod(bcellseq-1-1,3)+1;
rcellseq = mod(bcellseq-2-1,3)+1;
fstart = cut{seq(fcellseq)}(2);
bstart = cut{seq(bcellseq)}(2);
bstop = fstop
fstop = cut{seq(fcellseq)}(1);
cutr = 1;
rrecord = brecord;
% Clear the variable "fstoparray"
clear fstoparray;
[fstop, bstop, brecord, fstoparray, allarea] = genmixedcell(trial,fcellseq,bcellseq,rcellseq,fstart,fstop,bstart,bstop,cutr,rrecord,color,poly,nvpoly,cut,cutedge,seq,edges,allarea);
% Flag 9
% Swap the front cell and the remain cell and do mixed area computation again
temp = fcellseq;
fcellseq = rcellseq;
rcellseq = temp;
fstart = cut{seq(fcellseq)}(1);
rrecord = fstop;
rrecordarray = fstoparray;
allarea = special_genmixedcell(fcellseq,bcellseq,rcellseq,fstart,bstart,rrecord,rrecordarray,color,poly,nvpoly,cut,cutedge,seq,edges,allarea);

% Flag 10
bcellseq = 1;
fcellseq = mod(bcellseq-1-1,3)+1;
rcellseq = mod(bcellseq-2-1,3)+1;
fstart = cut{seq(fcellseq)}(1);
bstart = cut{seq(bcellseq)}(1);
bstop = fstop
%fstop = keepbstop
fstop = cut{seq(fcellseq)}(2);
cutr = 2;
rrecord = brecord;
% Clear the variable "fstoparray"
clear fstoparray;
[fstop, bstop, brecord, fstoparray, allarea] = genmixedcell(trial,fcellseq,bcellseq,rcellseq,fstart,fstop,bstart,bstop,cutr,rrecord,color,poly,nvpoly,cut,cutedge,seq,edges,allarea);
% Flag 11
% Swap the front cell and the remain cell and do mixed area computation again
temp = fcellseq;
fcellseq = rcellseq;
rcellseq = temp;
fstart = cut{seq(fcellseq)}(2);
rrecord = fstop;
rrecordarray = fstoparray;
allarea = special_genmixedcell(fcellseq,bcellseq,rcellseq,fstart,bstart,rrecord,rrecordarray,color,poly,nvpoly,cut,cutedge,seq,edges,allarea);

% Flag 12
bcellseq = 3;
fcellseq = mod(bcellseq-1-1,3)+1;
rcellseq = mod(bcellseq-2-1,3)+1;
fstart = cut{seq(fcellseq)}(2);
bstart = cut{seq(bcellseq)}(2);
bstop = fstop
fstop = cut{seq(fcellseq)}(1);
cutr = 1;
rrecord = brecord;
% Clear the variable "fstoparray"
clear fstoparray;
[fstop, bstop, brecord, fstoparray, allarea] = genmixedcell(trial,fcellseq,bcellseq,rcellseq,fstart,fstop,bstart,bstop,cutr,rrecord,color,poly,nvpoly,cut,cutedge,seq,edges,allarea);

mixedarea = dot([0 0 0 1 1 1],allarea);
clc;
title('Mixed Subdivision of Three Polygons P, Q and R');
xlabel(sprintf('Mixed area is %d + %d + %d = %d',allarea(4),allarea(5),allarea(6),mixedarea));
disp(sprintf('The mixed area is %d \n',allarea(4)+allarea(5)+allarea(6)));
disp(sprintf('Cell Area information:'));
disp(sprintf('======================'));
disp(sprintf('Blue: %d',allarea(1)));
disp(sprintf('Red: %d',allarea(2)));
disp(sprintf('Green: %d',allarea(3)));
disp(sprintf('Pink: %d',allarea(6)));
disp(sprintf('Cyan: %d',allarea(5)));
disp(sprintf('Yellow: %d',allarea(4)));
disp(sprintf('\n'));

matsize = size(lattice,2);
if (matsize~=0)
   % Initialize the graph of sparse resultant matrix
   % Plot on a new graph
   figure(2);
   % Set axis properties
   axis ij;
   axis square;
   axis([0 matsize 0 matsize]);
   hold on;
   axis off;
   
   sortlattice;
   resultant(P,Q,R);
   
   % Draw the boundary of the resultant matrix
   plot([0 matsize matsize 0 0], [0 0 matsize matsize 0], 'k');
   % Title of the graph
   title('Structure of the Sparse Resultant Matrix');
end

function [hull, numverts, edges, verts] = analysis(supp)

% Use "convhull" command to generate the convex hull.
% Note that it has some bugs in "convhull" commands, later, we will fix it.
hull = convhull(supp(1,:), supp(2,:));
% Since one connecting vertex is repeated, so we ignore it
numverts = size(hull, 2) - 1;
% Note that # vertices = # edges
% Store all the vector of edges of the convex hull 
% in an array "edges" (in anti-clockwise direction)
edges = zeros(2, numverts);
edges(:, 1:numverts) = supp(:, hull(2:numverts+1))...
   -supp(:, hull(1:numverts));
% Store all the vertices of the convex hull 
% in an array "verts" (in anti-clockwise direction)
verts = zeros(2, numverts);
verts(:, 1:numverts) = supp(:, hull(1:numverts));

function [index1a, index1b, index2a, index2b, index3a, index3b] = selectcut(verts1, verts2, verts3)

index1a = 0;
index1b = 0;
index2a = 0;
index2b = 0;
index3a = 0;
index3b = 0;
bigarea = 0;
numverts1 = size(verts1, 2);
numverts2 = size(verts2, 2);
numverts3 = size(verts3, 2);

for i = 1:numverts1-1
   for j = i+1:numverts1
      for k = 1:numverts2-1
         for l = k+1:numverts2
            for m = 1:numverts3-1
               for n = m+1:numverts3
                  vector1 = verts1(:,i) - verts1(:,j);
                  vector2 = verts2(:,k) - verts2(:,l);
                  vector3 = verts3(:,m) - verts3(:,n);
                  temparea = abs(det([vector1, vector2]))+abs(det([vector2, vector3]))+abs(det([vector3, vector1]));
                  if (temparea > bigarea)
                     bigarea = temparea;
                     index1a = i;
                     index1b = j;
                     index2a = k;
                     index2b = l;
                     index3a = m;
                     index3b = n;
                  end
               end
            end
         end
      end
   end
end

function d = genindex(a,b,c)
if (a>b)
   d = [mod([a-1:(b+c-1)], c)];
   d = d + ones(size(d, 1));
else
   d = [a:b];
end

function [fstop, bstop, brecord, fstoparray, allarea] = genmixedcell(t, f, b, r, fstart, fstop, bstart, bstop, cutr, rrecord, color, poly, nvpoly, cut, cutedge, seq, edges, allarea)
global mixed;
TRUE = 1;
FALSE = 0;

i = fstart
starte1 = edges{seq(f)}(:,i)
j = bstart
starte2 = -edges{seq(b)}(:,mod(j-1-1,nvpoly(seq(b)))+1)
overlap = FALSE;

while (det([starte2,starte1])>0)
   e1 = edges{seq(f)}(:,i)
   e2 = -edges{seq(b)}(:,mod(j-1-1,nvpoly(seq(b)))+1)
   pos = 0;
   while (det([e2,e1])>0)
      % FORM MIXED CELL
      tempmc{seq(f)} = [i mod(i+1-1,nvpoly(seq(f)))+1]
      tempmc{seq(b)} = [mod(j-1-1,nvpoly(seq(b)))+1 j]
      
      if ((i == fstop)|(j == bstop)|(overlap == TRUE))
         i
         fstop
         j
         bstop
         overlap
         fstoparray
         overlap = TRUE;
         tempmc{seq(r)} = rrecord
      else
         tempmc{seq(r)} = cut{seq(r)}(cutr)
      end
      
      if t==FALSE
         allarea = colorcell(tempmc, poly, seq, color, seq(r), seq(f), seq(b), allarea, 1);
         mixed = {mixed{:},{tempmc{:}}}
      end
         
      j = mod(j-1-1,nvpoly(seq(b)))+1;
      e2 = -edges{seq(b)}(:,mod(j-1-1,nvpoly(seq(b)))+1);
      brecord = j;
      pos = pos+1;
      fstoparray(pos) = mod(i+1-1,nvpoly(seq(f)))+1;
   end
      
   overlap = FALSE;
   
   if i == fstart
      bstoptemp = j;
   end
  
   i = mod(i+1-1,nvpoly(seq(f)))+1;
   starte1 = edges{seq(f)}(:,i);
   j = bstart;
   starte2 = -edges{seq(b)}(:,mod(j-1-1,nvpoly(seq(b)))+1);
end
fstop = i
fstoparray
bstop = bstoptemp

% Since the 2 unmixed cells is not adjacent, so the mixed area overlap on each other
function allarea = special_genmixedcell(f, b, r, fstart, bstart, rrecord, rrecordarray, color, poly, nvpoly, cut, cutedge, seq, edges, allarea)
global mixed;
i = fstart
starte1 = edges{seq(f)}(:,i)
j = bstart
starte2 = -edges{seq(b)}(:,mod(j-1-1,nvpoly(seq(b)))+1)

while (det([starte2,starte1])>0)
   e1 = edges{seq(f)}(:,i)
   e2 = -edges{seq(b)}(:,mod(j-1-1,nvpoly(seq(b)))+1)
   pos = 1
   while (det([e2,e1])>0)
      % FORM MIXED CELL
      tempmc{seq(f)} = [i mod(i+1-1,nvpoly(seq(f)))+1]
      tempmc{seq(b)} = [mod(j-1-1,nvpoly(seq(b)))+1 j]
      tempmc{seq(r)} = rrecordarray(pos);
      mixed = {mixed{:},{tempmc{:}}}
      pos = pos+1
      allarea = colorcell(tempmc, poly, seq, color, seq(r), seq(f), seq(b), allarea, 1);
      
      j = mod(j-1-1,nvpoly(seq(b)))+1
      e2 = -edges{seq(b)}(:,mod(j-1-1,nvpoly(seq(b)))+1)
   end
   i = mod(i+1-1,nvpoly(seq(f)))+1;
   starte1 = edges{seq(f)}(:,i)
   j = bstart
   starte2 = -edges{seq(b)}(:,mod(j-1-1,nvpoly(seq(b)))+1)
end

function chr = ccolor(n)
if n==1
   chr='b'
else
   if n==2
      chr='r'
   else
      chr='g'
   end
end

function allarea = colorcell(cell, poly, seq, color, n, nplus1, nplus2, allarea, celltype)
global eps;
global lattice_info;

xord = [poly{nplus1}(1,cell{nplus1}(1:2))+poly{nplus2}(1,cell{nplus2}(2))+poly{n}(1,cell{n}),...
      poly{nplus1}(1,cell{nplus1}(2:-1:1))+poly{nplus2}(1,cell{nplus2}(1))+poly{n}(1,cell{n})]
yord = [poly{nplus1}(2,cell{nplus1}(1:2))+poly{nplus2}(2,cell{nplus2}(2))+poly{n}(2,cell{n}),...
      poly{nplus1}(2,cell{nplus1}(2:-1:1))+poly{nplus2}(2,cell{nplus2}(1))+poly{n}(2,cell{n})]
fill(xord+eps(1), yord+eps(2), color{n+3});
allarea(n+3) = allarea(n+3) + polyarea(xord, yord)

intpt = polyintpt([xord+eps(1);yord+eps(2)],celltype,n+3)
if (intpt>0)
   formcell = zeros(1,3)
   formcell = [size(cell{1},2),size(cell{2},2),size(cell{3},2)]
   [u,v] = min(formcell)
   xord = poly{v}(1,cell{v})
   yord = poly{v}(2,cell{v})
   li = [xord, yord, v]'*ones(1,intpt)
   lattice_info = [lattice_info,li]
end

xord = [poly{nplus1}(1,cell{nplus1}(1:2))+poly{nplus2}(1,cell{nplus2}(1))+poly{n}(1,cell{n})]
yord = [poly{nplus1}(2,cell{nplus1}(1:2))+poly{nplus2}(2,cell{nplus2}(1))+poly{n}(2,cell{n})]
plot(xord+eps(1), yord+eps(2), ccolor(nplus1));
xord = [poly{nplus1}(1,cell{nplus1}(1:2))+poly{nplus2}(1,cell{nplus2}(2))+poly{n}(1,cell{n})]
yord = [poly{nplus1}(2,cell{nplus1}(1:2))+poly{nplus2}(2,cell{nplus2}(2))+poly{n}(2,cell{n})]
plot(xord+eps(1), yord+eps(2), ccolor(nplus1));
xord = [poly{nplus1}(1,cell{nplus1}(1))+poly{nplus2}(1,cell{nplus2}(1:2))+poly{n}(1,cell{n})]
yord = [poly{nplus1}(2,cell{nplus1}(1))+poly{nplus2}(2,cell{nplus2}(1:2))+poly{n}(2,cell{n})]
plot(xord+eps(1), yord+eps(2), ccolor(nplus2));
xord = [poly{nplus1}(1,cell{nplus1}(2))+poly{nplus2}(1,cell{nplus2}(1:2))+poly{n}(1,cell{n})]
yord = [poly{nplus1}(2,cell{nplus1}(2))+poly{nplus2}(2,cell{nplus2}(1:2))+poly{n}(2,cell{n})]
plot(xord+eps(1), yord+eps(2), ccolor(nplus2));

% celltype=0 means unmixed cell
% celltype=1 means mixed cell
function intpt = polyintpt(poly, celltype, cellinfo)
global nverts;
global mcptlst;
global umcptlst;
global cellptlst;
global lattice;

poly
size(poly,2)
if (size(poly,2)>2)
 index = convhull(poly(1,:), poly(2,:));
 poly = poly(:,index(:,1:size(index,2)-1));
 nverts = size(poly,2);

 [xmin,xminindex] = min(poly(1,:));
 [xmax,xmaxindex] = max(poly(1,:));

 % Note that we need to search from bottom edge --> top edge
 % We need to choose the x-min. and y-min element
 xplus1 = plus1(xminindex);
 while (poly(1,xplus1) == xmin)
    xminindex = xplus1;
    xplus1 = plus1(xplus1);
 end

 % We need to choose the x-max. and y-min element
 xminus1 = minus1(xmaxindex);
 while (poly(1,xminus1) == xmax)
    xmaxindex = xminus1;
    xminus1 = minus1(xminus1);
 end
 
 % Begin the big algorithm for finding the interior point in polygon
 botindex = xminindex;
 topindex = minus1(xminindex);
 botinv = [poly(:, botindex), poly(:, plus1(botindex))];
 topinv = [poly(:, plus1(topindex)), poly(:, topindex)];
 x = ceil(botinv(1,1));

 % Check whether the first top edge is a vertical line or not
 if (topinv(1,1)==topinv(1,2))
    % Change upper interval
    topindex = minus1(topindex);
    topinv = [poly(:, plus1(topindex)), poly(:, topindex)];
 end
   
 % Count the number of integer point in the polygon
 numofpts = 0;
 
 % Initialize intpt
 intpt = 0;
 
 % Initialize latpt
 latpt = [];
 
 while (1==1)
    botinv;
    topinv;
    if (x<=botinv(1,2))
       if (ceil(f(x,botinv))==floor(f(x,topinv)))
          % Store 1 point only
          y=ceil(f(x,botinv));
          plotdot(x, y, celltype, cellinfo);
          if (celltype == 0)
             umcptlst = [umcptlst, [x;y]];
          else
              mcptlst = [mcptlst, [x;y]];
          end
          numofpts = numofpts+1;
          cellptlst{cellinfo} = [cellptlst{cellinfo},[x;y]]
          lattice = [lattice,[x;y]]
       else
          % Store all the y points within the interval
          for y=ceil(f(x, botinv)):floor(f(x,topinv))
             plotdot(x, y, celltype, cellinfo);
             if (celltype == 0)
                umcptlst = [umcptlst, [x;y]];
             else
                mcptlst = [mcptlst, [x;y]];
             end
             numofpts = numofpts+1;
             cellptlst{cellinfo} = [cellptlst{cellinfo},[x;y]]
             lattice = [lattice,[x;y]]
          end
       end
       x=x+1;
    end
    if (x>botinv(1,2))
       % Change the lower interval
       botindex = plus1(botindex);
       botinv = [poly(:, botindex), poly(:, plus1(botindex))];
       if (botindex==xmaxindex)
          % Exit loop
          break
       end
    end
    if (x>topinv(1,2))
       % Change upper interval
       topindex = minus1(topindex);
       topinv = [poly(:, plus1(topindex)), poly(:, topindex)];
       if (plus1(topindex)==xmaxindex)
          % Exit loop
          break
       end
    end
 end
 intpt = numofpts
end
      
function plotdot(x,y,i,t)
sym = '';
brg = ['b', 'r', 'g'];
x,y,i,t
% If the lattice point in unmixed cell, plot in 'x' sign.
% If the lattice point in mixed cell, plot in '+' sign.
if (i==0)
   sym = 'kx';
else
   sym = strcat(brg(t-3),'+');
end
plot(x,y,sym);

function v1 = plus1(v0)
global nverts;
v1 = mod(v0+1-1, nverts)+1;

function v1 = minus1(v0)
global nverts;
v1 = mod(v0-1-1, nverts)+1;

function y=f(x,inv)
diff = inv(:,2)-inv(:,1);
m = diff(2)/diff(1);
y = m*x+(inv(2,1)-m*inv(1,1));

function resultant(P,Q,R)
global lattice;
global lattice_info;

fid = fopen('res.m','w');

fprintf(fid,'%% Lattice points are:\n');

for i = 1:size(lattice, 2)
   fprintf(fid,'%% (%d,%d)\n',lattice(:,i));
end

fprintf(fid,'\n');
fprintf(fid,'%% Enter your coefficient here (Only input coeff. of f1, f2)\n');
fprintf(fid,'%% Format: coeff{1} = [1st 2nd 3nd];\n');
fprintf(fid,'\n');
fprintf(fid,'%% f3 is used for helping to solve f1=f2=0\n');
fprintf(fid,'%% Note the first coefficient of f3 is ignored\n');
fprintf(fid,'\n');

matsize = (size(lattice,2));
% Compute M = M0 + uM1, where M0 and M1 are all coefficients of the polynomial systems
fprintf(fid,'M0 = zeros(%d);\n',matsize);
fprintf(fid,'M1 = zeros(%d);\n',matsize);
for i = 1:matsize
   latpt = lattice(:,i);
   formpt = lattice_info(1:2,i);
   f = lattice_info(3,i);
   if (f==1)
      polygon = P;
   elseif (f==2)
      polygon = Q;
   else
      polygon = R;
   end
   new_polygon = [latpt-formpt]*ones(1,size(polygon,2))+polygon;
   for t = 1:size(polygon,2)
      j = find((lattice(1,:)==new_polygon(1,t))&(lattice(2,:)==new_polygon(2,t)));
      % Note that the first coefficient of f3 must be u, where u is a parameter. 
      if (f==3)
         if (t==1)
            fprintf(fid,'M1(%d,%d) = 1;\n',i,j);
         else
            fprintf(fid,'M0(%d,%d) = coeff{%d}(%d);\n',i,j,f,t);
         end
      else
         fprintf(fid,'M0(%d,%d) = coeff{%d}(%d);\n',i,j,f,t);
      end
      matgraph(i,j,f);
   end
end

fprintf(fid,'\n');
fprintf(fid,'%% Solve the generalized Eigenvalue problem.\n');
fprintf(fid,'%% Sample code:\n');
fprintf(fid,'%% [v, d] = eig(M0, M1);\n');
fprintf(fid,'\n');
fprintf(fid,'%% Use the degree in lattice points to find x and y...\n');
fprintf(fid,'%% Sample code:\n');
fprintf(fid,'%% x = [];\n');
fprintf(fid,'%% y = [];\n');
fprintf(fid,'%%\n');
fprintf(fid,'%% for u = 1:size(M0, 1)\n');
fprintf(fid,'%%    x = [x; v(1,u)/v(2,u)]\n');
fprintf(fid,'%%    y = [y; v(3,u)/v(4,u)]\n');
fprintf(fid,'%% end\n');
fprintf(fid,'\n');
fclose(fid);

function sortlattice
global lattice;
global lattice_info;

% Use the built-in function sortrows
% Initially, our format of lattice is [x0 x1 ... ]
%                                     [y0 y1 ... ],
% The format of lattice need to transpose in order it can use sortrows.

% Ordering of each pair (x,y): Sort by y-cord in ascending order,  
%    if y-cord is the same, then sort by x-cord in ascending order.
[useless, sortindex] = sortrows(lattice([2 1],:)');

new_lattice = lattice(:,sortindex');
new_lattice_info = lattice_info(:,sortindex');

lattice = new_lattice;
lattice_info = new_lattice_info;

function matgraph(i, j, c)
color{1} = 'b';
color{2} = 'r';
color{3} = 'g';

% Note i here stands for row and j stands for column.
% However, in cordinate, i will become y where j will become x
fill([j-1, j, j, j-1], [i-1, i-1, i, i], color{c});
plot([j-1, j], [i-1, i-1], color{c});
plot([j, j], [i-1, i], color{c});
plot([j-1, j], [i, i], color{c});
plot([j-1, j-1], [i-1, i], color{c});
