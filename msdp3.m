% (8/11/99)
%
%    This Matlab code computes a mixed
% subdivision of the convex hulls of THREE input
% supports (using the ``maxpinch'' heuristic),
% and then illustrates the output with color
% graphics.  The heuristic is described in the fourth
% paragraph following this one.  Also, the
% program gives some information, roughly comparing
% the sizes of three different resultants.
%
% NOTE: Numerical output is mostly given as an
%       INTEGER.  This is not a problem for inputs
%       with INTEGER coordinates.  If
%       floating point output is desired,
%       then one can simply change
%       the %d's (in the sprintf's)
%       to, say, %15.15g
%
%    The program does a little input cleaning,
% computes the convex hulls (correcting a
% bug in Matlab 5.0's convhull command in the
% process), and then computes the mixed
% subdivision.
%
%   Computing the mixed subdivision entails
% walking around the second polygon, until a
% vertex v2 is found which can be properly
% attached to the ``first'' vertex of
% the first polygon.  (Proper attachment
% means that the corresponding angle cones
% intersect only at their vertices.)
% This portion takes worst-case time
% O(m).
%
%    The heuristic which continues the mixed
% subdivision computation consists of a walk
% around both polygons, going through all
% proper attachments, until the
% pinch area is maximized.  (The pinch
% area of two polygons touching at a
% single vertex is simply the sum of
% the area(s) of the parallelogram(s)
% defined by the complements of the
% two angle cones.  So the 4 (or 2) edges
% form half the edges of the ``pinched''
% parallelograms.)  This walk around
% the polygons takes an extra O(m+n)
% time, but will hopefully help later
% in constructing sparse resultant
% matrices with simpler structure.
%
%    The subdivision computation is then completed
% by ``expanding a shock wave'' emanating
% from the second polygon, constantly keeping
% track of the subdivision edges which
% lie in the convex sum.  The loops
% (toward the end) attempt to perform
% this expansion efficiently, so
% as to obtain an O(mn) worst-case complexity
% (after the two polygons have been
% properly placed) with small constant (around
% 1/4 on average).
%
% NOTE:  The subdivision thus computed
%        is found deterministically.
%        The only thing special about
%        the subdivisions computed with
%        this code is that they seek
%        to maximize the sum of the
%        areas of two particular cells,
%        in the (rough) hopes of minimizing
%        the final number of mixed cells.
%
%    The output consists of the mixed area
% of the two polygons (new=nan indicates
% an error), any necessary error messages,
% and a drawing of the mixed subdivision.
% In particular, one usually sees polygon
% with minimal number of vertices in red
% the other in blue, and all the
% mixed cells in pink.  (In the degenerate
% instances where one of the polygons is an
% edge or a point, some of the polygons might
% not be drawn.)  The union of
% all the colored polygons is precisely
% the Minkowski sum of the two polygons.
 
 
% suppi should be an input matrix of size 2 by
% n_i, containing the x and y coordinates of
% the ith support set of interest.
%
% NOTE1: I should really do a little extra
%        optimizing for the case of
%        exactly two edges, when those
%        edges are affinely independent
%        (a better version of folding)...
%
% NOTE2: This code does not do lattice optimization,
%        i.e., if every thing is actually living on
%        a sublattice of finite index, this will
%        *NOT* be taken into account.
%

function [areas,mixedcells]=msdp3(supp1,supp2,supp3)
       
% Let's find if there are wrong dimensions, repeats, or 
% singletons in our first or second supports...
supp1=check(supp1,1); supp2=check(supp2,2); 
supp3=check(supp3,3);  

% Set up the Newton polygons for our supports... 
[numverts(1) newt{1} frontedge{1}]=polygon(supp1); 
[numverts(2) newt{2} frontedge{2}]=polygon(supp2); 
[numverts(3) newt{3} frontedge{3}]=polygon(supp3); 

for i=1:3
 polys(i)=0; 
 if numverts(i)>2 
  polys(i)=polyarea(newt{i}(1,:),newt{i}(2,:)); 
 end
end

% Compute the classical and multihomogeneous Bezout 
% bounds...
[bez,multi,bmink,mmink]=bezout(newt{1},newt{2},newt{3}); 

color{1}='b'; color{2}='r'; color{3}='g'; % Polygon 
                                          % colors...

% Cell colors...
ccolor{1}=[1   1   .6 ]; 
ccolor{2}=[.8  1   1  ]; 
ccolor{3}=[1   .9  1  ];  

bolor(1)=3;   bolor(2)=1;   bolor(3)=2;  % Indices for 
                                         % cell colors...

% What's missing from the pair {1,2,3}\{i}?  
% (This also helps keep track of the swaps...)
pair(1)=1;    pair(2)=2;    pair(3)=3; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort the three polygons according to # of vertices...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numverts(2)<numverts(3) 
 swap=pair(2);  pair(2)=pair(3);   pair(3)=swap; 
end
if numverts(pair(1))<numverts(pair(2)) 
 swap=pair(2);  pair(2)=pair(1);   pair(1)=swap; 
end
if numverts(pair(2))<numverts(pair(3)) 
 swap=pair(2);  pair(2)=pair(3);   pair(3)=swap; 
end

special=2; 
% If any hull is empty, or any size is wrong, stop!
if numverts(pair(3))==0 
 % Later on, this WILL include some plotting... !!!
 mixedcells={{[],[]},{[],[]},{[],[]}}; 
 lengths=zeros(1,3);
 if numverts(pair(2))==0
  if numverts(pair(1))~=1 
   disp(sprintf('Toric Resultant Vanishes Identically!'));
   special=-1; 
   areas(1)=nan; areas(2)=nan; areas(3)=nan; 
   lengths(1)=nan; lengths(2)=nan; lengths(3)=nan; 
   lengtht=nan; 
  else 
   special=0; 
   areas(1)=0; areas(2)=0; areas(3)=0; 
   lengths(pair(1))=1; lengtht=1; 
   disp(sprintf('Toric Resultant Degree is %d+%d+%d = %d\n\n',...
         lengths(1),lengths(2),lengths(3),lengtht));
  end  
 elseif numverts(pair(2))==1
  areas(1)=0; areas(2)=0; areas(3)=0; 
  if numverts(pair(1))==1
   disp(sprintf('Toric Resultant Degree is 0')); 
   special=0;  
   lengtht=0;  
  else 
   lengths(pair(2))=1; lengtht=1;  
   disp(sprintf('Toric Resultant Degree is %d+%d+%d = %d\n\n',...
         lengths(1),lengths(2),lengths(3),lengtht));
  end
 elseif numverts(pair(2))==2
  if numverts(pair(1))==2
   if det([frontedge{pair(1)}(:,1),frontedge{pair(2)}(:,1)])==0
    ed1=frontedge{pair(1)}(:,1); 
    ed2=frontedge{pair(2)}(:,1); 
    % Let's find the normalization for the 
    % lattice generated by the two vectors 
    % ed1 and ed2...
    g1=abs(gcd(ed1(1),ed1(2))); 
    g2=abs(gcd(ed2(1),ed2(2))); 
    lengths(pair(1))=g2; 
    lengths(pair(2))=g1; 
    lengtht=dot(lengths,[1 1 1]); 
    disp(sprintf('Toric Resultant Degree is %d+%d+%d = %d\n',...
         lengths(1),lengths(2),lengths(3),lengtht));
    special=1; 
   else 
    disp(sprintf('Toric Resultant Vanishes Identically!'));
    special=-1;
    areas(1)=nan; areas(2)=nan; areas(3)=nan;
    lengths(1)=nan; lengths(2)=nan; lengths(3)=nan;
    lengtht=nan;
   end 
  else 
   disp(sprintf('Toric Resultant Vanishes Identically!'));
   special=-1;
   areas(1)=nan; areas(2)=nan; areas(3)=nan;
   lengths(1)=nan; lengths(2)=nan; lengths(3)=nan;
   lengtht=nan; 
  end
 else 
  disp(sprintf('Toric Resultant Vanishes Identically!'));
  special=-1;
  areas(1)=nan; areas(2)=nan; areas(3)=nan;
  lengths(1)=nan; lengths(2)=nan; lengths(3)=nan;
  lengtht=nan; 
 end
 return 
end

% If, by perverse chance, we have at least one 
% one point support...
if numverts(pair(3))==1
 lengths=zeros(1,3);
 if numverts(pair(2))==1
  special=0;
  lengths(1)=0; lengths(2)=0; lengths(3)=0;
  lengtht=0;
  if numverts(pair(1))==1
   disp(sprintf('There''s nothing to draw!\n')); 
  end
 else 
  special=0; 
  lengths(pair(3))=1; lengtht=1; 
  if all([(numverts(pair(2))==2),(numverts(pair(1))==2),...
          (det([frontedge{pair(1)}(:,1),...
                frontedge{pair(2)}(:,2)])==0)])
     lengths(pair(3))=0; lengtht=0; 
  end
 end
end

% Check if at least one of the polygons 
% is an edge (this is simpler since we've just sorted)... 
doubled=0; 
if numverts(pair(3))==2 % The only worrisome cases are 
                      % 432, 322, and 222...
 if numverts(pair(2))>2 
  % This means we have just one edge, 
  % so let's do a little optimization...
  rot=[0 1;-1 0]*frontedge{pair(3)}(:,1);   
  dots1=zeros(1,numverts(pair(1))); 
  for i=1:numverts(pair(1))
   dots1(i)=dot(newt{pair(1)}(:,i),rot); 
  end
  dots2=zeros(1,numverts(pair(2))); 
  for i=1:numverts(pair(2))
   dots2(i)=dot(newt{pair(2)}(:,i),rot); 
  end
  height1=max(dots1)-min(dots1);  
  height2=max(dots2)-min(dots2);  
  if height1>height2 
   swap=pair(2);  pair(2)=pair(1);   pair(1)=swap;
  end
 else

  % Try to move earlier parallel pairs to the second 
  % and third positions...
  if numverts(pair(1))==2
   if det([frontedge{pair(1)}(:,1),frontedge{pair(2)}(:,1)])==0
    swap=pair(3);  pair(3)=pair(1);   pair(1)=swap;
   end
   if det([frontedge{pair(1)}(:,1),frontedge{pair(3)}(:,1)])==0
    swap=pair(2);  pair(2)=pair(1);   pair(1)=swap;
   end 
  end

  % Now that any REALLY degenerate pairs are at the 
  % end, let's check if there are any...
  if det([frontedge{pair(2)}(:,1),frontedge{pair(3)}(:,1)])==0
   % This means newt{pair(2)} and newt{pair(3)} are parallel 
   % edges... 

   % Set up the edge colors for the first polygon...
   fillcolor1=color{pair(1)};
   color{pair(1)}=zeros(1,numverts(pair(1)));
   color{pair(1)}(:)=fillcolor1;

   % No need to worry about these hulls, since  
   % we're no longer treating a 3 by 2 resultant...
   % (But this should be changed later when 
   %  we include the Sylvester resultant into 
   %  all of this...)
   hulls{pair(2)}=[]; 
   hulls{pair(3)}=[]; 
   hull23={[],[]};     % No need to draw the boundary of 
                       % the sum...
   mixedcells{pair(1)}={[],[]}; % Likewise, there are no {2,3}-mixed cells...
   [foo,vert2]=cut(newt{pair(2)},frontedge{pair(2)},2,...
                   newt{pair(3)},frontedge{pair(3)});
   hift=newt{pair(2)}(:,vert2)-newt{pair(3)}(:,1);
   shift=zeros(2,2); 
   shift(1,:)=hift(1); shift(2,:)=hift(2);
   newt23=[newt{pair(2)}(:,3-vert2),(newt{pair(3)}(:,2)+hift)]; 
   front23=newt23(:,2)-newt23(:,1); 
   [start,finish]=cut(newt{pair(1)},frontedge{pair(1)},...
                      numverts(pair(1)),newt23,front23); 
   
   doubled=1; % This means we've encountered 
              % one of the worst degeneracies... 

   if and(numverts(pair(1))==2,...
          det([frontedge{pair(1)}(:,1),...
               frontedge{pair(2)}(:,1)])==0) 
    % If we have three parallel edges...
    % (and we have a higher order resultant...)
    hulls{pair(1)}=[]; % No first hull to draw...
    hull123={[newt{pair(1)}(:,[start,finish]),...
       (newt{pair(1)}(:,finish)+frontedge{pair(3)}(:,1)),...
       (newt{pair(1)}(:,finish)+front23(:)),...
       (newt{pair(1)}(:,finish)+frontedge{pair(3)}(:,1)),...
       newt{pair(1)}(:,finish)],...
       [color{pair(1)}(1),color{pair(3)},color{pair(2)},...
        color{pair(2)},color{pair(3)},color{pair(1)}(1)]};
    mixedcells={{[],[]},{[],[]},{[],[]}}; 
    % Set the *NORMALIZED* toric resultant parameters! 
    special=0; 
    lengths=zeros(1,3); 
    lengtht=0; 
   else 
    % This means newt{pair(2)} and newt{pair(3)} are 
    % parallel edges but newt{pair(1)} is NOT parallel to 
    % newt{pair(2)} or newt{pair(3)}...

    % Are these endpoints adjacent?... 
    slide3=frontedge{pair(3)}(:,1); 
    if finish==(mod(start,numverts(pair(1)))+1)
     % Use the first polygon with its color...
     reorder=mod([1:numverts(pair(1))]+finish-2,...
                 numverts(pair(1)))+1;
     hulls{pair(1)}=newt{pair(1)}; 
     % A dirty fix for an orientation bug (7/30/99)...
     % (The case where there is a single 12 cell 
     %  and a single 13 cell, attached to the 
     %  end of newt{pair(1)}, in such a way that the 
     %  ``ccw'' is ruined...)    
     if det([front23,...
             (newt{pair(1)}(:,finish)-...
              newt{pair(1)}(:,start))])<0
      slide3=-slide3; 
      front23=-front23; 
     end 
     hull123={[newt{pair(1)}(:,reorder(:)),...
      (newt{pair(1)}(:,reorder(end))+slide3),...
      (newt{pair(1)}(:,reorder(end))+front23),...
      (newt{pair(1)}(:,reorder(1))+front23),...
      (newt{pair(1)}(:,reorder(1))+slide3)],...
      [color{pair(1)}(1:(end-1)),color{pair(3)},...
       color{pair(2)},color{pair(1)}(1),...
       color{pair(2)},color{pair(3)}]}; 
    else
     % Set the indices for the two pieces of
     % the cut polygon...
     if finish>start
      good=[start:finish];
      bad=mod([(finish-1):(numverts(pair(1))+start-1)],...
              numverts(pair(1)));
      bad=bad+1;
     else
      good=mod([(start-1):(numverts(pair(1))+finish-1)],...
               numverts(pair(1)));
      good=good+1;
      bad=[finish:start];
     end
     shift=zeros(2,size(good,2));
     shift(1,:)=front23(1); 
     shift(2,:)=front23(2);
     hulls{pair(1)}={(newt{pair(1)}(:,good)+shift),...
                     newt{pair(1)}(:,bad)};
     hull123={[(newt{pair(1)}(:,good)+shift(:,:)),...
                (newt{pair(1)}(:,good(end))+...
                  slide3),...
                 newt{pair(1)}(:,bad),...
                 (newt{pair(1)}(:,bad(end))+...
                  slide3)],...
               [color{pair(1)}(good(1:(end-1))),...
                color{pair(2)},color{pair(3)},...
                color{pair(1)}(bad(1:(end-1))),...
                color{pair(3)},color{pair(2)}]};
    end
    % Now define the cells that were cut into (or appended
    % to) the first polygon...
    paralleloa=zeros(2,4); parallelob=zeros(2,4); 
    paralleloa(1,:)=newt{pair(1)}(1,start);
    paralleloa(2,:)=newt{pair(1)}(2,start);
    paralleloa(:,2)=paralleloa(:,2)+...
                    newt{pair(1)}(:,finish)-...
                    newt{pair(1)}(:,start);
    paralleloa(:,3:4)=[paralleloa(:,2),paralleloa(:,4)]+...
                      [slide3,slide3];
    parallelob(1,:)=newt{pair(1)}(1,start)+slide3(1);  
    parallelob(2,:)=newt{pair(1)}(2,start)+slide3(2);  
    parallelob(:,2)=parallelob(:,2)+...
                    newt{pair(1)}(:,finish)-...
                    newt{pair(1)}(:,start);  
    parallelob(:,3:4)=[parallelob(:,2),parallelob(:,4)]+...
                      [(front23-slide3),(front23-slide3)];
    mixedcells{pair(2)}={[paralleloa],...
                    [color{pair(1)}(1),color{pair(3)}]};
    mixedcells{pair(3)}={[parallelob],...
                    [color{pair(1)}(1),color{pair(2)}]}; 

    % Set the *NORMALIZED* toric resultant parameters!
    special=1; 
    lengths=zeros(1,3);
    ed2=frontedge{pair(2)}(:,1);
    ed3=frontedge{pair(3)}(:,1);
    % Let's find the normalization for the
    % lattice generated by the two vectors
    % ed2 and ed3...
    g2=abs(gcd(ed2(1),ed2(2)));
    g3=abs(gcd(ed3(1),ed3(2)));
    lengths(pair(2))=g3;
    lengths(pair(3))=g2;
    lengtht=dot(lengths,[1 1 1]);
   end
  % NOTE: If newt{pair(2)} and newt{pair(3)} are edges which 
  %       are NOT parallel, nothing need be done... 
  else 
   % This means we have at least two edges,  
   % and if there are three, they are 
   % pairwise affinely independent...
   if numverts(pair(1))==2 
    % This is the latter case...  

    % Set up the edges ed2 and ed3, 
    % so that they go in ccw order... 
    ed2=frontedge{pair(2)}(:,1); 
    ed3=frontedge{pair(3)}(:,1); 

    % Ensure that ed2 to ed3 is ccw...
    if det([ed2,ed3])<0
     swap=ed2; ed2=ed3; ed3=swap; 
     swap=pair(2);  pair(2)=pair(3);   pair(3)=swap; 
    end 

    % Now let's throw in the final edge...
    ed1=frontedge{pair(1)}(:,1); 
    % Flip ed1 if lies in the cone generated by 
    % ed2 and ed3...
    if and(det([ed2,ed1])>0,det([ed1,ed3])>0)
     ed1=-ed1; 
    end

    % Now we must eliminate the possibility 
    % ed1, ed2, ed3 lying in a common half-space...
    % (The first determinant ensures the cone 
    %  is strongly convex.)
    if all([det([ed2,ed1])>0,...
            det([ed2,ed3])>0,det([ed3,ed1])>0])
     ed3=-ed3; % Flip ed3 if it's between ed2 and ed1...
     % ...and correct the order as well! 
     swap=ed1; ed1=ed3; ed3=swap; 
     swap=pair(1);  pair(1)=pair(3);   pair(3)=swap; 
    end
    % (The first determinant ensures the cone 
    %  is strongly convex.)
    if all([det([ed1,ed3])>0,...
            det([ed1,ed2])>0,det([ed2,ed3])>0])
     ed2=-ed2; % Flip ed2 if it's between ed1 and ed3...
     % ...and correct the order as well! 
     swap=ed1; ed1=ed2; ed2=swap; 
     swap=pair(1);  pair(1)=pair(2);   pair(2)=swap; 
    end

    % Finally, we can define the 23 cell...
    para23=zeros(2,4);
    para23(:,2)=ed2; 
    para23(:,3)=ed2+ed3; 
    para23(:,4)=ed3;
    hull23={para23,...
            [color{pair(2)},color{pair(3)},...
             color{pair(2)},color{pair(3)}]}; 
    mixedcells{pair(1)}={[para23],[color{pair(2)},color{pair(3)}]};

    % ...and the 31 and 12 cells...
    para31=zeros(2,4); 
    para31(:,2)=ed3; 
    para31(:,3)=ed3+ed1; 
    para31(:,4)=ed1;
    para12=zeros(2,4); 
    para12(:,2)=ed1; 
    para12(:,3)=ed1+ed2; 
    para12(:,4)=ed2;

    % The next few lines were changed (8/8/99) 
    % to fix a bug in delta(.) (the edges *ARE* 
    % needed in order to figure out the rays...)
    hulls{pair(2)}=[[0;0],ed2];  
    hulls{pair(3)}=[[0;0],ed3]; 
    hulls{pair(1)}=[[0;0],ed1]; 

    hull123={[ed2,(ed2+ed3),ed3,(ed3+ed1),...
              ed1,(ed1+ed2)],...
             [color{pair(3)},color{pair(2)},color{pair(1)},...
              color{pair(3)},color{pair(2)},color{pair(1)}]};
    mixedcells{pair(2)}={[para31],[color{pair(3)},color{pair(1)}]};
    mixedcells{pair(3)}={[para12],[color{pair(1)},color{pair(2)}]}; 
    doubled=1;  
   end
  end
 end % Every other case can be left alone!... 
end

if not(doubled) % If not already handled by an  
% exceptional case of doubled edges, 
% compute the Minkowski sum of convex hulls of the 
% first two (colored) supports, the mixed area, and then 
% plot a mixed subdivision derived via the 
% max-pinch heuristic...
% NOTE:  This function implicitly assumes that the 
%        sets of edge colors of the two input polygons 
%        are disjoint.
 [hulls23,hull23,mixcells23]=...
        mink(newt{pair(2)},frontedge{pair(2)},color{pair(2)},...
             newt{pair(3)},frontedge{pair(3)},color{pair(3)},...
             '');

 % This is the most interesting part:  Compute the edges 
 % of the first Minkowski sum, and then combine it with 
 % the third polygon!...
 frontedge23=[hull23{1}(:,2:end),hull23{1}(:,1)]-...
              hull23{1}(:,:); 

 if size(hull23{1},2)<=2 % Put any potentially degenerate 
                         % polygon second...  
  [hulls123,hull123,mixcells123]=...
      mink(newt{pair(1)},frontedge{pair(1)},color{pair(1)},...
           hull23{1},frontedge23,hull23{2},...
           color{pair(2)}(1)); 
  % Define the output cells at last...
  mixedcells{pair(1)}={[],[]}; 
  mixedcells{pair(2)}={[],[]};
  mixedcells{pair(3)}=mixcells123{2}; 

  % ...and the degenerate ones too... 
  hulls{pair(3)}=[]; 
  hulls{pair(2)}=hulls123{2}; 
  hulls{pair(1)}=hulls123{1};
 else 
  [hulls123,hull123,mixcells123]=...
      mink(hull23{1},frontedge23,hull23{2},...
           newt{pair(1)},frontedge{pair(1)},color{pair(1)},...
           color{pair(2)}); 
  % Define the output cells at last...
  mixedcells{pair(1)}=mixcells23{1};
  mixedcells{pair(2)}=mixcells123{1}; 
  mixedcells{pair(3)}=mixcells123{2}; 

  % ...and the degenerate ones too... 
  for i=1:2
   hulls{pair(i+1)}=hulls23{i};  
  end
  % Find the nonempty summand of (newt2+newt3)+newt1...
  thenon=hulls123{2}; 
  if isempty(thenon) 
   thenon=hulls123{1}; 
  end
  hulls{pair(1)}=thenon;   
 end 
end

% subplot(1,1,1); % This gets UNcommented when used for real 
                  % (not for cards) 
newplot; hold on; 

% Fill (if allowed) the colored Newton polygons...
for i=1:3
 if size(hulls{i},1)==1 
  fill(hulls{i}{1}(1,:),hulls{i}{1}(2,:),color{i}(1));
  fill(hulls{i}{2}(1,:),hulls{i}{2}(2,:),color{i}(1));
 else 
  if size(hulls{i},2)>2
   fill(hulls{i}(1,:),hulls{i}(2,:),color{i}(1));
  end
 end
end

% Draw the boundary of the ``first'' Minkowski sum hull!...
if size(hull23{1},2)>2
 drawboundary(hull23);  
end

% Draw the boundary of the final Minkowksi sum...
drawboundary(hull123);  

% Draw the mixed cells and compute the mixed 
% cell areas, type by type...
for i=1:3
 areas(i)=drawcells(mixedcells{i},ccolor{i}); 
end

areat=dot(areas,[1 1 1]); 
minky=areat+dot(polys,[1 1 1]); 
if special==2 
 disp(sprintf('Toric Resultant Degree is %d+%d+%d = %d, Minkowski Sum Area = %.1f\n\n',...
              areas(1),areas(2),areas(3),areat,minky));
 xlabel(sprintf('Toric Resultant Degree is %d+%d+%d = %d, Minkowski Sum Area = %.1f\n(Classical=%d, Multihomogeneous=%d, CMink=%.1f, MMink=%d)',...
              areas(1),areas(2),areas(3),areat,...
              minky,bez,multi,bmink,mmink));

elseif special==-1 
 disp(sprintf('Toric Resultant Vanishes Identically!'));
else 
 disp(sprintf('Toric Resultant Degree is %d+%d+%d = %d, Minkowski Sum Area = %.1f\n\n',...
         lengths(1),lengths(2),lengths(3),lengtht,minky));
 xlabel(sprintf('Toric Resultant Degree is %d+%d+%d = %d, Minkowski Sum Area = %.1f\n(Classical=%d, Multihomogeneous=%d, CMink=%.1f, MMink=%d)',...
              lengths(1),lengths(2),lengths(3),lengtht,...
              minky,bez,multi,bmink,mmink));
end

hold off; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function below, which computes the 
% *colored* Minkowski sum of two ordered 
% simple polygons, is the heart of 
% new(a,b) 
% NOTE1:  The function assumes newt1 and newt2 
%        are, respectively, 2 by k1 and 2 by k2 
%        arrays, WITH ANY DEGENERATE POLYGONS 
%        OCCURING SECOND (!!) AND ANY POLYGON 
%        WITH 2 COLORS FIRST (!!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [hulls,minkhull,mixcells]=...
         mink(newt1,frontedge1,color1,... 
              newt2,frontedge2,color2,goesfirst)

numverts1=size(newt1,2);
numverts2=size(newt2,2);

% If the ``goesfirst'' color code is empty, 
% make it something nonempty, but still 
% nonsensical...
if isempty(goesfirst)
 goesfirst='x';
end

% If the color code has only one character, 
% make all edge colors equal...
if size(color1,2)==1 
 fillcolor1=color1;  
 color1=zeros(1,numverts1); 
 color1(:)=fillcolor1;
end

% newt2 is assumed to be monochromatic, so 
% we only need do the following...
fillcolor2=color2;          
color2=zeros(1,numverts2); 
color2(:)=fillcolor2;

% Now that we have the Newton polygons (with their color 
% codes) is either one a point?  Draw 
% something and stop if either one is...
if numverts2==1
 hulls{2}=[];  
 minkhull={newt1,color1}; 
 mixcells={{[],[]},{[],[]}};
 if numverts1==1
  hulls{1}=[]; 
 else 
  hulls{1}=newt1; 
 end
 return
end 

% If either Newton polygon is an edge, then this also  
% simplifies matters...
if numverts2==2 % (the second polygon is an edge)
 % Get ready to set up the cells...
 shift=zeros(2,2);  
 second=frontedge2(:,1); 
 hulls{2}=newt2; % Another fix (8/8/99) 
                 % to make sure the 
                 % rays in delta(.) are 
                 % well-defined...
 
 % ...but this depends on the disposition of the edge(s).
 if numverts1==2 % The first polygon is also an edge...
  first=frontedge1(:,1);
  hulls{1}=newt1; % Another fix (8/8/99) 
                  % to make sure the 
                  % rays in delta(.) are 
                  % well-defined...
  area=abs(det([first,second])); 
  if area==0 % They're both parallel...
   [start,finish]=cut(newt1,first,2,newt2,second); 
   hift=newt1(:,finish)-newt2(:,1);
   shift=[hift,hift];
   minkhull={[newt1(:,[start,finish]),...
              (newt2(:,[2,1])+shift)],...
	     [color1(1),color2(1),color2(2),color1(2)]}; 
   mixcells={{[],[]},{[],[]}};
  else % They sum up to form a nondegenerate cell...
   parallelo=zeros(2,4);
   parallelo(1,1:2)=newt1(1,:); 
   parallelo(2,1:2)=newt1(2,:);
   parallelo(:,3:4)=[parallelo(:,2),parallelo(:,1)]+...
                    [frontedge2(:,1),frontedge2(:,1)]; 
   % Let's ensure that this new hull is oriented 
   % counter-clockwise...
   if det([frontedge1(:,1),frontedge2(:,1)])<0 
    parallelo=parallelo(:,end:-1:1); 
   end 
   minkhull={parallelo,...
             [color1(1),color2(1),color1(2),color2(2)]}; 
   mixcells={{[],[]},{[],[]}};
   if or(goesfirst==color1(1),goesfirst==color2(1))
    mixcells{2}={[parallelo],[color1(1),color2(1)]};
   else 
    mixcells{1}={[parallelo],[color1(1),color2(1)]};
   end
  end
  return % Ah...  We're returning from mink into msdp3...
 else % First polygon has at least 3 vertices (which 
      % usually implies positive area), but the second is 
      % just an edge... 

  [start,finish]=cut(newt1,frontedge1,numverts1,...
                     newt2,frontedge2); 

  slide2=frontedge2(:,1); 
  % Are these endpoints adjacent?...
  if finish==(mod(start,numverts1)+1)
   % Use the first polygon with its color...
   reorder=mod([1:numverts1]+finish-2,numverts1)+1; 
   hulls{1}=newt1; 
   % A dirty fix for an orientation bug (7/30/99)...
   if det([slide2,...
           (newt1(:,finish)-newt1(:,start))])<0
    slide2=-slide2; 
   end 
   minkhull={[newt1(:,reorder(:)),...
	     (newt1(:,reorder(end))+slide2),...
	     (newt1(:,reorder(1))+slide2)],...
	     [color1(1:(end-1)),color2(1),color1(1),...
              color2(1)]};
  else
   % Set the indices for the two pieces of
   % the cut polygon...
   if finish>start
    good=[start:finish];
    bad=mod([(finish-1):(numverts1+start-1)],numverts1);
    bad=bad+1;
   else
    good=mod([(start-1):(numverts1+finish-1)],numverts1);
    good=good+1;
    bad=[finish:start];
   end
   shift=zeros(2,size(good,2));
   shift(1,:)=slide2(1); shift(2,:)=slide2(2);
   hulls{1}={(newt1(:,good)+shift),newt1(:,bad)}; 
   minkhull={[(newt1(:,good)+shift),newt1(:,bad)],...
	     [color1(good(1:(end-1))),color2(1),...
			    color1(bad(1:(end-1))),...
                            color2(2)]};
  end
  % Now define the cell that was cut into (or appended 
  % to) the first polygon...
  parallelo=zeros(2,4);
  parallelo(1,:)=newt1(1,start);
  parallelo(2,:)=newt1(2,start);
  parallelo(:,2)=parallelo(:,2)+newt1(:,finish)-...
                                newt1(:,start);
  parallelo(:,3:4)=[parallelo(:,2),parallelo(:,4)]+...
                   [slide2,slide2];
  mixcells={{[],[]},{[],[]}}; 
  if or(goesfirst==color1(1),goesfirst==color2(1))
   mixcells{2}={[parallelo],[color1(1),color2(1)]}; 
  else 
   mixcells{1}={[parallelo],[color1(1),color2(1)]}; 
  end
  return
 end
else % The hardest part is when each Newton polygon 
     % has at least 3 vertices (which usually 
     % implies each has positive area)... 
     % (NOTE:  Although this *should* work for 
     %         the case of newt1 being two adjacent
     %         parallel edges of different colors, 
     %         it does NOT.  So I remove this 
     %         exception BEFORE I call this 
     %         function (7/29/99)...) 

 % Set up a good initial vert(1)
 v1=1; v1b=numverts1; 
 % Hopefully, the 4 lines below fix a  
 % bug I noticed... (parallel edges in 
 % Minkowski sums causing problems?) 
 done=0; % Part of a new bug fix (7/28/99): 
         % What if newt1 is a SERIES of 
         % parallel edges?... 
 while and(det(frontedge1(:,[v1b,v1]))==0,not(done)) 
  v1=mod(v1,numverts1)+1; 
  v1b=mod(v1b,numverts1)+1; 
  if v1==1
   done=1; 
  end 
 end

 % Get ready to search for an initial vert(2)
 v2=1; v2b=numverts2; 
 % Hopefully, the 4 lines below fix a  
 % bug I noticed... (parallel edges in 
 % Minkowski sums causing problems?) 
 while det(frontedge2(:,[v2b,v2]))==0 
  v2=mod(v2,numverts2)+1; 
  v2b=mod(v2b,numverts2)+1; 
 end

 % Let's find the orientations of the Newton polygons... 
 ccw=[sign(det(frontedge1(:,[v1b,v1]))) ,...  
      sign(det(frontedge2(:,[v2b,v2])))];  
 if ccw(1)==0
  ccw(1)=1; 
 end
 % The above pathology is assumed NOT to 
 % occur in newt2... (7/28/99)

 % Now let's find a v2 on newt2 whose angle 
 % cone is disjoint from the angle cone of 
 % the vertex 1 on newt1...  
 % (The main test is to see if a frontedgei(vert(i)+1)  
 % (or -frontedgei(vert(i)) of polygon i is within the 
 % cone generated by frontedgej(vert(j)+1) and 
 % -frontedgej(vert(j)).) There are four ``AND'' 
 % conditions to be ``OR'ed''...
 % NOTE: Sign of polygon j determines the factor in the inequality...
 while conesect(frontedge1(:,[v1b v1]),ccw(1),...
                frontedge2(:,[v2b v2]),ccw(2)) 
    v2=mod(v2,numverts2)+1 ;  
    v2b=mod(v2b,numverts2)+1 ; 
 end  
end

vert=[v1 v2]; 

% Now that we have a single proper placement, let's 
% run through *all* of them... 

% Let's set up orientations for the indexing of our 
% pinched cells...
% (Depending on the orientations of the Newton polygons, 
%  we'll want to compute max(det(x,y),0)+max(det(u,v),0), 
%  where [x y] = [ front1 back2 ] or [front1 front2] 
%                [ back1  back2 ] or [back1  front2], 
%        according as the signs are [1 1], [1 -1], [-1 1], and [-1 -1]
%  and similarly for 
%        [u v] = [back1 front2], [back1 back2], [front1 front2], 
%                or [front1 back2] . )
move2=ccw(1)*ccw(2) ; % The way to walk on the second Newton polygon...

% Let's get ready...
area=-1;  

paralleloa=zeros(2,4); parallelob=zeros(2,4);
done=0; 
% Let's try to find a vertex to vertex placement of 
% the Newton polygons, so that the pinched cell  
% areas are maximized...
newv2=vert(2); % Another fix on 7/29/99... 
               % (This one actually stopped 
               %  the unmixed case of three 
               %  unit squares from failing!)
while not(done) 
 worked=0; cycled=0; 
 % NB1: The while loop below has the opposite condition 
 %      as our initial search for a good initial position...
 % NB2: Note that we now move both v1 and v2... 
 while not(... 
        or(...
         conesect(frontedge1(:,[v1b,v1]),ccw(1),...
                  frontedge2(:,[v2b,v2]),ccw(2)),...
         cycled)) % The last bit takes care of a pathology I 
                  % blundered over...
                  % (If the second polygon is skinny, and 
                  %  things are just right (or wrong), it's 
                  %  possible to cycle around the second 
                  %  polygon forever!)

  % Let's set up our pinched cells...
  [pincha,pinchb]=pinch(frontedge1(:,[v1b,v1]),ccw(1),...
                        frontedge2(:,[v2b,v2]),ccw(2)); 
  
  % Sum the areas of the pinched cells...
  areaa=det(pincha); areab=det(pinchb);
  newarea=max(areaa,0)+max(areab,0); 

  % Are these pinched cells better?... 
  if or(newarea>area,and(newarea==area,or(areaa<=0,areab<=0)))
   opt1=v1; opt2=v2; 
   area=newarea; 
  end

  % If this is the first run through the inner loop,
  % keep track of this v2...
  if worked==0
   firstv2=v2; 
  end

  % This v2 worked, so let's note this:
  worked=1;  

  % Update the vertex on the second Newton polygon...
  v2=mod(v2+move2-1,numverts2)+1; 
  v2b=mod(v2-2,numverts2)+1;
  % Hopefully, the 4 lines below fix a  
  % bug I noticed... (parallel edges in 
  % Minkowski sums causing problems?) 
  while det(frontedge2(:,[v2b v2]))==0 
   v2=mod(v2+move2-1,numverts2)+1; 
   v2b=mod(v2-2,numverts2)+1; 
  end
   
  % ...and make sure it hasn't cycled around!
  if v2==firstv2
   cycled=1;
  end

  % Check if we've already wrapped around... 
  if and(v1==vert(1),v2==vert(2))
   done=1; 
  end
 end

 % Update the vertices... 
 if and(worked,not(done))  
  v1=mod(v1,numverts1)+1; 
  v1b=mod(v1-2,numverts1)+1; 

  v2=mod(v2-move2-1,numverts2)+1;
  v2b=mod(v2-2,numverts2)+1;  
  % Hopefully, the 4 lines below fix a 
  % bug I noticed... (parallel edges in 
  % Minkowski sums causing problems?) 
  while det(frontedge2(:,[v2b v2]))==0 
   v2=mod(v2-move2-1,numverts2)+1; 
   v2b=mod(v2-2,numverts2)+1; 
  end
  newv2=v2; 
 else 
  v2=mod(v2+move2-1,numverts2)+1 ;  
  v2b=mod(v2-2,numverts2)+1;  
  % Hopefully, the 4 lines below fix a 
  % bug I noticed... (parallel edges in 
  % Minkowski sums causing problems?) 
  while det(frontedge2(:,[v2b v2]))==0 
   v2=mod(v2+move2-1,numverts2)+1; 
   v2b=mod(v2-2,numverts2)+1; 
  end
  if newv2==v2 
   v1=mod(v1,numverts1)+1;
   v1b=mod(v1-2,numverts1)+1;
  end
 end

 % Check again if we've already wrapped around... 
 if and(v1==vert(1) , v2==vert(2))
  done=1;  
 end

end 

% Let's prepare for our output arrays... 
shift=zeros(2,numverts2) ;
hift=newt1(:,opt1)-newt2(:,opt2) ;
shift(1,:)=hift(1); shift(2,:)=hift(2);
opt1b=mod(opt1-2,numverts1)+1; opt2b=mod(opt2-2,numverts2)+1;

% Let's set up the two Newton polygons and their 
% fill colors...
hulls{1}=newt1; 
hulls{2}=(newt2+shift);  

mixcells={{[],[]},{[],[]}}; % ...and fill the output arrays.
minkhull={[],[]}; %
parallelo=zeros(2,4); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's finally define the convex hull of the Minkowski
%% sum...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE:  It appears that I'm implictly assuming that the 
%        the first polygon (when both are positive dimensional)
%        is oriented COUNTERCLOCKWISE.  This OK, modulo 
%        a bug-fix I just made for the case of two 
%        linearly independent edges (above, 7/22/99). 
move=zeros(2,2); % move(``half'' , polygon ) 
if move2==1             % Depending on the orientations,
 edgeind2=[opt2b opt2]; % the indices of the starting edges of 
 move=[1 -1 ; -1 1];    % the ``forward'' and ``backward'' 
else                    % parallelograms will be different...
 edgeind2=[opt2 opt2b]; % (We'll also need to know which way the 
 move=[1 1 ; -1 -1];    % indices move...)
end                       
edgeind1=[opt1 opt1b];  

% Let's first find the edges of the 
% Minkowski sum by ``merging'' the edges 
% of our two Newton polygons.  This entails 
% a geometric sort...
limit={[],[]}; 
part={[],[]};
pcolor={[],[]}; 
for i=1:2
 % Let's first find the forward (or backward)
 % limit of the edges on the second
 % Newton polygon, from the pinch...
 lim=0;
 initedge1=frontedge1(:,edgeind1(i));
 e2=edgeind2(i); 
 while (ccw(1)*move(i,2)*det([frontedge2(:,e2),initedge1]))>0;
  e2=mod(e2+move(i,2)-1,numverts2)+1;
  lim=lim+1;
 end
 if lim>0
  lim=lim-1; % The limit was too big by one...
  % Define the first two endpoints of 
  % part i of the convex hull of the Minkowski sum...
  limit{i}=[limit{i} lim];
  hift=newt1(:,opt1)-newt2(:,opt2);  
  startind2=mod(opt2+(lim+1)*move(i,2)-1,numverts2)+1;
  start2=newt2(:,startind2);
  qart=[start2 (start2+ccw(1)*move(i,1)*initedge1)]+...
       [hift hift]; 
  part{i}=[part{i} qart]; 
  pcolor{i}=[pcolor{i} color1(edgeind1(i))]; 
   
  % Now that we have our first limit, let's start
  % the sort... 
  e1=mod(edgeind1(i)+move(i,1)-1,numverts1)+1;
  edge1=frontedge1(:,e1); 
  initedge2=frontedge2(:,edgeind2(i));
  finaledge2=frontedge2(:,mod(edgeind2(i)+...
                              lim*move(i,2)-1,...
                              numverts2)+1);
  while (ccw(1)*move(i,2)*det([initedge2,edge1]))>0
   % Do we need to lower the limit?...
   while (ccw(1)*move(i,2)*det([finaledge2,edge1]))<=0
    lim=lim-1;
    finaledge2=frontedge2(:,mod(edgeind2(i)+...
                                lim*move(i,2)-1,...
                                numverts2)+1); 
   end
   startind2=mod(opt2+(limit{i}(end)+1)*move(i,2)-1,...
                 numverts2)+1;
   hift=part{i}(:,end)-newt2(:,startind2); 
   for j=limit{i}(end):-1:(lim+1) 
    pcolor{i}=[pcolor{i} color2(startind2)];  
    startind2=mod(startind2-move(i,2)-1,numverts2)+1;
    part{i}=[part{i},(newt2(:,startind2)+hift)]; 
   end
   part{i}=[part{i},(part{i}(:,end)+ccw(1)*move(i,1)*edge1)];
   pcolor{i}=[pcolor{i},color1(e1)];  
   limit{i}=[limit{i} lim];

   e1=mod(e1+move(i,1)-1,numverts1)+1;
   edge1=frontedge1(:,e1);
  end
  startind2=mod(opt2+(limit{i}(end)+1)*move(i,2)-1,...
                numverts2)+1;
  hift=part{i}(:,end)-newt2(:,startind2); 
  for j=limit{i}(end):-1:0 
   pcolor{i}=[pcolor{i},color2(startind2)];
   startind2=mod(startind2-move(i,2)-1,numverts2)+1;
   part{i}=[part{i},(newt2(:,startind2)+hift)]; 
  end 
 end
end 

% Starting from the pinch, incrementally store the 
% parallelograms, first on one side, and then on the other 
% side, of the first Newton polygon.
for i=1:2 
 e1=edgeind1(i); v1=opt1; 
 edge1=ccw(1)*move(i,1)*frontedge1(:,e1); 
 initedge2=move(i,2)*frontedge2(:,edgeind2(i));
 thelimit=limit{i}; 
 % Since we already have our limits, let's start 
 % incrementally drawing our parallelograms...

 for j=1:size(thelimit,2)
  e2=edgeind2(i); edge2=initedge2; 
  parallelo(1,:)=newt1(1,v1);
  parallelo(2,:)=newt1(2,v1);
  for k=0:thelimit(j)
   parallelo(:,2)=parallelo(:,1)+edge1; 
   parallelo(:,3:4)=[parallelo(:,2),parallelo(:,4)]+...
                    [edge2,edge2]; 
   if or(goesfirst==color1(e1),goesfirst==color2(e2))
    mixcells{2}{1}=[mixcells{2}{1},parallelo];
    mixcells{2}{2}=[mixcells{2}{2},[color1(e1),color2(e2)]];
   else 
    mixcells{1}{1}=[mixcells{1}{1},parallelo];
    mixcells{1}{2}=[mixcells{1}{2},[color1(e1),color2(e2)]];
   end
   parallelo(:,1)=parallelo(:,1)+edge2; 
   e2=mod(e2+move(i,2)-1,numverts2)+1; 
   edge2=move(i,2)*frontedge2(:,e2); 
  end
  v1=mod(v1+move(i,1)-1,numverts1)+1;  
  e1=mod(e1+move(i,1)-1,numverts1)+1;
  edge1=ccw(1)*move(i,1)*frontedge1(:,e1);
 end
end

% Define the length (and indices) of the portion of the second 
% polygon lying on the Minkowksi hull...
lim1=-1;                      % NOTE:  In this case (both newt1 and newt2
if not(isempty(limit{1}))     %        having positive area), it is 
 lim1=limit{1}(1);            %        impossible for both limit{1} and 
end                           %        limit{2} to be empty... 
lim2=-1;
if not(isempty(limit{2}))
 lim2=limit{2}(1); 
end
long2=numverts2-lim1-lim2-3; 
if long2==-1 % Delete the start of part{2} if necessary...
 part{2}=part{2}(:,2:end);
end

part{2}=part{2}(:,end:-1:1);     % Turn around part{2} and pcolor{2}
pcolor{2}=pcolor{2}(end:-1:1);   % because of the way our ``shockwave''
			         % from newt2 was constructed...

hift=newt1(:,opt1)-newt2(:,opt2); % The second polygon is shifted a bit...

firstvert2=lim2+1; 
if isempty(part{1}) 
 long2=long2+1;            % If part{1} is empty, we'll need an extra vertex...
end                     

if isempty(part{2})
 firstvert2=-lim1-long2-2; % If instead, part{2} is empty, then we'll 
end                        % index a little differently... 
 

% Let's FINALLY get the indexing correct...
ind2=mod(opt2+move(2,2)*([1:long2]+firstvert2)-1,numverts2)+1;
ind2b=mod(opt2+move(2,2)*([0:long2]+firstvert2)-1,numverts2)+1;
shift=zeros(2,size(ind2,2));
shift(1,:)=hift(1); shift(2,:)=hift(2); 
part2=newt2(:,ind2)+shift;  
if isempty(part{1})
 pcolor2=color2(ind2b(1:(end-1))); % pcolor2 will need one color less 
else                               % iff part{1} is empty... 
 pcolor2=color2(ind2b);            % (Since we're walking about CCW...) 
end

% Define the length (and indices) of the portion of the first polygon 
% lying on the Minkowksi hull...
long1=numverts1-size(limit{1},2)-size(limit{2},2)-1;  
if long1==-1 % Delete the end of part{1} and pcolor{1} if necessary...
 part{1}=part{1}(:,1:(end-1));
end

firstvert1=size(limit{1},2);
if isempty(part{2}) 
 long1=long1+1;            % If part{2} is empty, we'll need an extra vertex...
end                     

if isempty(part{1})
 firstvert1=-size(limit{2},2)-long1-1; % If instead, part{1} is empty, 
end                                    % then we'll index a bit differently... 

ind1=mod(opt1+([1:long1]+firstvert1)-1,numverts1)+1;
ind1b=mod(opt1+([0:long1]+firstvert1)-1,numverts1)+1;
part1=newt1(:,ind1);
if isempty(part{2})
 pcolor1=color1(ind1b(1:(end-1))); % pcolor1 will need one color less 
else                               % iff part{2} is empty... 
 pcolor1=color1(ind1b);            % (Since we're walking about CCW...) 
end

% Define the colored edge list of the Minkowski sum 
% at last!...
if isempty(part{1})
 minkhull{1}=[part{2},part2,part{1},part1]; 
 minkhull{2}=[pcolor{2},pcolor2,pcolor{1},pcolor1]; 
else 
 minkhull{1}=[part{1},part1,part{2},part2];    
 minkhull{2}=[pcolor{1},pcolor1,pcolor{2},pcolor2]; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here are some functions I use in new...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This one's for error-handling in the input...
function newsupp=check(supp,ind);  
s=size(supp,2) ; 

if size(supp,1)~=2 
 disp(sprintf(['\nYour ', card(ind),...
               ' support has the wrong size!\n\n']))
 newsupp=[];
 return
end

if s==1 
 newsupp=supp; 
else 
 newsupp=unique(supp','rows')'; 
end

% This one's for cardinalities in English...
function word=card(ind);
if ind==1
 word='first';
elseif ind==2
 word='second';
elseif ind==3
 word='third'; 
else 
 word='nth?'; 
end

% This one's for setting up the Newton polygon 
% (including the edge list) of a support... 
function [numverts,newt,frontedge]=polygon(supp) 
s=size(supp,2); 
if s<=1 
 numverts=s; frontedge=[]; 
 if s==0 
  newt=[];
 else 
  newt=supp;
 end
 return 
elseif s==2
 newt=supp; numverts=2; 
 frontedge=supp(:,2)-supp(:,1); 
 frontedge=[frontedge -frontedge];  
else 
 % Set up the vectors emanating from the first point...
 bump=zeros(2,s-1); 
 bump(1,1:s-1)=supp(1,1); 
 bump(2,1:s-1)=supp(2,1); 
 bumped=zeros(2,s-1) ; 
 bumped(:,1:s-1)=supp(:,2:s)-bump(:,1:s-1) ; 
 dets=zeros(1,s-2) ; 
 for i=1:s-2
  dets(i)=(det([bumped(:,1),bumped(:,i+1)])==0) ;  
 end
 % Check if the convex hull is just an edge...
 if all(dets) 
  numverts=2; true=[1 2]; 
  dmin=dot(bumped(:,1),supp(:,1)) ; dmax=dmin ; 
  amin=1 ; amax=1 ; 
  % Find the extreme points in the direction of the 
  % edge to define the edge...
  for i=2:s
   if dot(bumped(:,1),supp(:,i))>dmax  
    dmax=dot(bumped(:,1),supp(:,i)); 
    amax=i ; 
   end
   if dot(bumped(:,1),supp(:,i))<dmin  
    dmin=dot(bumped(:,1),supp(:,i)); 
    amin=i ; 
   end
  end
  newt=[supp(:,amin),supp(:,amax)];
  frontedge=supp(:,amax)-supp(:,amin); 
  frontedge=[frontedge -frontedge]; 
 else % If the Newton polygon has positive area...
  % ...let's finally use the convex hull command.
  gull=convhull(supp(1,:),supp(2,:)) ;  
  pumverts=size(gull,2)-1 ; 

  % Let's (allocate and) define an array for the 
  % edges of the Newton polygon...
  grontedge=zeros(2,pumverts) ; 
  grontedge(:,1:pumverts)=supp(:,gull(2:pumverts+1))... 
                          -supp(:,gull(1:pumverts)) ; 
  % ...and let's check for repeats... 
  true=[];  
  for i=1:pumverts 
   if det([grontedge(:,mod(i-2,pumverts)+1), ...
           grontedge(:,mod(i-1,pumverts)+1)])~=0 
    true=[true,mod(i-1,pumverts)+1] ; 
   end
  end 
  % ...and modify as necessary.
  numverts=size(true,2);  
  if numverts<pumverts 
   % Throw out points on edge interiors...
   hull=zeros(1,numverts); 
   hull=gull(true(:));
   newt=[supp(1,hull(:));supp(2,hull(:))];  
   % Compress the edgelist...
   frontedge=zeros(2,numverts); 
   frontedge(1,:)=supp(1,hull([2:numverts,1]))-...
                   supp(1,hull(1:numverts));  
   frontedge(2,:)=supp(2,hull([2:numverts,1]))-...
                   supp(2,hull(1:numverts));  
  else 
   frontedge=grontedge; 
   newt=[supp(1,gull(1:pumverts));supp(2,gull(1:pumverts))];  
  end
 end
end 

% This is to see if an edge vector lies in the  
% angle cone defined by two consecutive edges...
% (Is inedge of polygon i within the cone generated by 
% -bedge and fedge?)  
% NOTE: The orientation of the polygon that bedge and 
%       fedge come from determines the sign factor in the 
%       inequality...
function yes=invector(inedge,bfedge,ccw)
 yes=and((ccw*det([bfedge(:,2),inedge]))>=0,...
         (ccw*det([bfedge(:,1),inedge]))>=0); 

% This one's for checking if the two angle 
% cones (defined by two pairs of consecutive 
% edges) intersect... 
% (This comes down to ``OR'ing'' four invector  
% conditions, but one of them can be ommitted 
% since there can never be exactly one true 
% condition.)  
function yes=conesect(bfedge1,ccw1,bfedge2,ccw2)
 yes=any([invector(bfedge2(:,2),bfedge1,ccw1),...
     invector(-bfedge2(:,1),bfedge1,ccw1),...
     invector( bfedge1(:,2),bfedge2,ccw2) ]);  

% This one's for building the pinched 
% parallelograms...
function [pincha,pinchb]=pinch(bfedge1,ccw1,bfedge2,ccw2)
 if ccw1*ccw2==1
  pincha=[-ccw2*bfedge2(:,1),ccw1*bfedge1(:,2)];
  pinchb=[-ccw1*bfedge1(:,1),ccw2*bfedge2(:,2)];
 else 
  pincha=[-ccw2*bfedge2(:,2),ccw1*bfedge1(:,2)];
  pinchb=[-ccw1*bfedge1(:,1),ccw2*bfedge2(:,1)];
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing the Bezout bounds...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bez,multi,bmink,mmink]=bezout(newt1,newt2,newt3) 
 disp(sprintf('Up to monomial factors, the (total degree,x-degree,y-degree) triples are...'));
 if isempty(newt1)
  disp(sprintf('Undefined!'));
  td1=nan; md1=[nan,nan];  
 else 
  % Find the degrees, shifts for cornering, the 
  % multidegrees, the Bezout numbers, and 
  % the multihomogeneous Bezout numbers...
  degs1=newt1(1,:)+newt1(2,:); % Find the total degrees
  minx1=min(newt1(1,:)); miny1=min(newt1(2,:));
  maxx1=max(newt1(1,:)); maxy1=max(newt1(2,:));
  td1=max(degs1)-minx1-miny1; 
  md1=[maxx1-minx1,maxy1-miny1];
  disp(sprintf('(%d,%d,%d)',td1,md1(1),md1(2)));
 end
 if isempty(newt2)
  disp(sprintf('Undefined!'));
  td2=nan; md2=[nan,nan]; 
 else 
  % Find the degrees, shifts for cornering, the 
  % multidegrees, the Bezout numbers, and 
  % the multihomogeneous Bezout numbers...
  degs2=newt2(1,:)+newt2(2,:); % of the monomial terms...
  minx2=min(newt2(1,:)); miny2=min(newt2(2,:));
  maxx2=max(newt2(1,:)); maxy2=max(newt2(2,:));
  td2=max(degs2)-minx2-miny2;
  md2=[maxx2-minx2,maxy2-miny2];
  disp(sprintf('(%d,%d,%d)',td2,md2(1),md2(2)));
 end
 if isempty(newt3)
  disp(sprintf('Undefined!\n'));
  td3=nan; md3=[nan,nan]; 
 else 
  % Find the degrees, shifts for cornering, the 
  % multidegrees, the Bezout numbers, and 
  % the multihomogeneous Bezout numbers...
  degs3=newt3(1,:)+newt3(2,:); % 
  minx3=min(newt3(1,:)); miny3=min(newt3(2,:));
  maxx3=max(newt3(1,:)); maxy3=max(newt3(2,:));
  td3=max(degs3)-minx3-miny3;
  md3=[maxx3-minx3,maxy3-miny3];
  disp(sprintf('(%d,%d,%d)\n',td3,md3(1),md3(2))); 
 end

 bezout1=td2*td3; 
 bezout2=td1*td3;
 bezout3=td1*td2;
 bez=bezout1+bezout2+bezout3;
 bmink=bez+(td1^2+td2^2+td3^2)/2;
 disp(sprintf('Dense Resultant Degree is %d+%d+%d = %d, Minkowski Sum Area = %.1f\n',...
 	      bezout1,bezout2,bezout3,bez,bmink));

 % Compute the multihomogeneous Bezout bound...
 mb3=md1(1)*md2(2)+md1(2)*md2(1);
 mb2=md1(1)*md3(2)+md1(2)*md3(1);
 mb1=md2(1)*md3(2)+md2(2)*md3(1);
 multi=mb3+mb2+mb1; 
 mmink=multi+md1(1)*md1(2)+md2(1)*md2(2)+md3(1)*md3(2);
 disp(sprintf(...
              'Multihomogeneous Resultant Degree is %d+%d+%d = %d, Matrix Size = %d\n',mb1,mb2,mb3,...
              multi,mmink));

% This one's for finding mixed subdivision of a 
% polygon and an edge, using the ``fold'' heuristic: 
% That is, one finds a good place to ``fold'' 
% (or ``cut'') the first polygon in two, 
% and then one inserts a mixed cell between 
% the two halve...
function [start,finish]=cut(newt1,frontedge1,numverts1,...
                            newt2,frontedge2) 
 % Find the orientation of the first polygon and 
 % then find the Minkowski sum in this simple 
 % special case... 
 v1b=numverts1; v1=1; 

 if numverts1==2  % A simplification (and bug fix) for 
                  % when newt1 is just an edge. (7/30/99)
  start=1; finish=2; 
  if and(det([frontedge1(:,1),frontedge2(:,1)])==0,...
         dot(frontedge1(:,1),frontedge2(:,1))<0) 
   % Things need only be changed if newt1 and newt2 
   % are parallel and point in the opposite direction...
   start=2; finish=1; 
  end 
  return 
 end

 % Hopefully, the 4 lines below fix a
 % bug I noticed... (parallel edges in
 % Minkowski sums causing problems?)
 done=0; % Part of a new bug fix (7/28/99):
         % What if newt1 is a SERIES of
         % parallel edges?...
 while and(det(frontedge1(:,[v1b,v1]))==0,not(done))
  v1=mod(v1,numverts1)+1;
  v1b=mod(v1b,numverts1)+1;
  if v1==1
   done=1;
  end
 end

 ccw(1)=sign(det(frontedge1(:,[v1b,v1])));  
 if ccw(1)==0 
  ccw(1)=1; 
 end % This will be thrown out soon...
     % Since the newt1=edge case will 
     % soon be ruled out by outer work... (7/29/99)
 
 % Find the first alternation between
 % frontedge2(:,1) being in or out...
 % (There was a bug here (7/11/99), but 
 % I fixed it by flipping frontedge2(:,1) 
 % if one can succesfully cycle around 
 % with no alternations.  This should 
 % optimized in a better way later...)
 last=invector(frontedge2(:,1),...
               frontedge1(:,[v1b,v1]),ccw(1));
 v1b=mod(v1b,numverts1)+1; v1=mod(v1,numverts1)+1; vi=v1; 
 repeat=0; 
 next=invector(frontedge2(:,1),...
               frontedge1(:,[v1b,v1]),ccw(1));
 while and(last==next,not(repeat))
  v1b=v1; v1=mod(v1,numverts1)+1;
  if v1==vi 
   repeat=1;
  end
  next=invector(frontedge2(:,1),...
                frontedge1(:,[v1b,v1]),ccw(1));
 end
 first=v1; 

 % (Here's the remainder of the bug fix.)
 if repeat
  frontedge2=-frontedge2;
  last=invector(frontedge2(:,1),...
                frontedge1(:,[v1b v1]),ccw(1));
  v1b=mod(v1b,numverts1)+1; v1=mod(v1,numverts1)+1; vi=v1;
  next=invector(frontedge2(:,1),...
                frontedge1(:,[v1b,v1]),ccw(1)); 
  while last==next 
   v1b=v1; v1=mod(v1,numverts1)+1;
   next=invector(frontedge2(:,1),...
                 frontedge1(:,[v1b v1]),ccw(1));
  end
  first=v1; 
 end

 % Find the second such alternation...
 last=next;
 v1b=v1; v1=mod(v1,numverts1)+1;
 next=invector(frontedge2(:,1),...
               frontedge1(:,[v1b v1]),ccw(1));
 while last==next
  v1b=v1; v1=mod(v1,numverts1)+1;
  next=invector(frontedge2(:,1),...
                frontedge1(:,[v1b,v1]),ccw(1));
 end
 second=v1;

 % Now set the ordered end points of the path
 % which can have frontedge2(:,1) emanating from
 % each vertex...
 if last % If first gave a vector IN the angle cone...
  start=second; finish=mod(first-2,numverts1)+1;
 else
  start=first; finish=mod(second-2,numverts1)+1;
 end

% Draw the boundary of the Minkowski hull!...
function drawboundary(hull); 
lim=size(hull{1},2);
if lim>1 
 for i=1:lim 
  plot(hull{1}(1,mod([i-1,i],size(hull{1},2))+1),...
       hull{1}(2,mod([i-1,i],size(hull{1},2))+1),...
       hull{2}(i));
 end
end 

% Draw the mixed cells (and compute the sum of 
% their areas)! 
function area=drawcells(mixcells,cellcolor)
area=0;
for j=0:(size(mixcells{1},2)/4-1)
 i=1+4*j;
 l=1+2*j; 
 leg1=mixcells{1}(:,i+1)-mixcells{1}(:,i);
 leg2=mixcells{1}(:,i+2)-mixcells{1}(:,i);
 area=area+abs(det([leg1,leg2]));
 fill(mixcells{1}(1,[i:(i+3)]),...
      mixcells{1}(2,[i:(i+3)]),...
      cellcolor);
 
 plot(mixcells{1}(1,[i,(i+1)]),...
      mixcells{1}(2,[i,(i+1)]),...
      mixcells{2}(l));
 plot(mixcells{1}(1,[(i+1),(i+2)]),...
      mixcells{1}(2,[(i+1),(i+2)]),...
      mixcells{2}(l+1));
 plot(mixcells{1}(1,[(i+2),(i+3)]),...
      mixcells{1}(2,[(i+2),(i+3)]),...
      mixcells{2}(l));
 plot(mixcells{1}(1,[(i+3),i]),...
      mixcells{1}(2,[(i+3),i]),...
      mixcells{2}(l+1));
end

