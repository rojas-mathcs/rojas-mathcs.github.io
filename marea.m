%    Simple Matlab code for computing the mixed area 
% of the convex hulls of a pair of input point sets. 
%
%    The program does a little input cleaning,
% computes the convex hulls (correcting a
% bug in Matlab 5.0's convhull command in the
% process), and then computes the mixed
% area.
%
%    The input checking takes O(m^2+n^2)
% worst-case time, where m (resp. n) is the
% minimum (resp. maximum) of the numbers of
% vertices in the convex hulls of the input
% point sets.  The quadratic complexity (as
% opposed to quasi-linear) is due to quickly
% written code (in the function check(...) ),
% and can certainly be improved.
%
%    Computing the mixed area entails
% walking around the polygon with fewer vertices, 
% until a vertex v2 is found which can be properly
% attached to the ``first'' vertex of
% the polygon with more vertices.  (Proper attachment
% means that the corresponding angle cones
% intersect only at their vertices.)
% This portion takes worst-case time
% O(m).
%
%    The mixed area computation then continues 
% with the computation of a mixed subdivision. 
% This is done by ``expanding a shock wave'' emanating
% from one of the polygons, constantly keeping
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
%        is found deterministically, 
%        but with no attempt at any 
%        sort of optimization.
%
%    The mixed area is then simply the 
% sum of the mixed cell areas, which 
% is computed along with the mixed subdivision. 
%
%    The output consists of the mixed area
% of the two polygons (msd=nan indicates
% an error), and any necessary error messages. 

% suppi is an input matrix of size 2 by numvertsi+1, 
% containing the x and y coordinates of 
% the ith support set of interest.
%
% NB:  The trickiest parts of writing this code, in  
%      ascending order, were 
%      (a) keeping track of where to put ... 
%          and + so matlab wouldn't surreptitiously 
%          do something strange. 
%      (b) correctly coding all the ``conical'' 
%          intersection conditions, 
%      (c) handling all degenerate instances, 
%      (d) indexing all the arrays correctly, and
%      (e) not making every conceivable dumb 
%          mistake possible!  : )
function area=marea(supp1,supp2)

% Let's find if there are wrong dimensions, repeats, or 
% singletons in our first or second supports...
supp1=check(supp1,1); supp2=check(supp2,2);  

% Set up the Newton polygons for our supports... 
[numverts1 frontedge1 hull1]=polygon(supp1); 
[numverts2 frontedge2 hull2]=polygon(supp2); 

% If either hull is empty, or of the wrong size, stop!
if or(numverts1==0,numverts2==0) 
 area=nan; 
 return
end

% Make the smaller one second...
if numverts2>numverts1 
 swap=numverts2; numverts2=numverts1; numverts1=swap;
 swap=hull2; hull2=hull1; hull1=swap;
 swap=frontedge2; frontedge2=frontedge1; frontedge1=swap; 
 swap=supp2; supp2=supp1; supp1=swap; 
end

% Now that we have the Newton polygons (with the one with 
% fewest vertices last), is either one a point?  Draw something 
% and stop if either one is...
if or(numverts1==1,numverts2==1)
 area=0; 
 return
end 

% If either Newton polygon is an edge, then this also  
% simplifies matters...
if numverts2==2 % (the second polygon is an edge)
 second=frontedge2(:,1); 

 if numverts1==2 % The first polygon is also an edge...
  first=frontedge1(:,1) ; 
  area=abs(det([first,second]));
  return
 else % First polygon has positive area, but the second is just an edge...
  % Find the orientation of the first polygon and 
  % then find the Minkowski sum in this simple 
  % special case...
  ccw(1)=sign(det(frontedge1(:,1:2))); 

  % Find the first alternation between
  % frontedge2(:,1) being in or out...
  % (There was a bug here (7/11/97), but
  % I fixed it by flipping frontedge2(:,1)
  % if one can succesfully cycle around
  % with no alternations.  This should
  % optimized in a better way later...)
  last=invector(frontedge2(:,1),frontedge1(:,[numverts1 1]),ccw(1));
  v1=2; repeat=0;
  next=invector(frontedge2(:,1),frontedge1(:,[1 2]),ccw(1));
  while and(last==next,not(repeat))
   v1b=v1; v1=mod(v1,numverts1)+1;
   if v1==2 
    repeat=1; 
   end
   next=invector(frontedge2(:,1),...
                 frontedge1(:,[v1b v1]),ccw(1));
  end
  first=v1;

  % (Here's the remainder of the bug fix.)
  if repeat
   frontedge2=-frontedge2;
   last=invector(frontedge2(:,1),...
                 frontedge1(:,[numverts1 1]),ccw(1));
   v1=2;
   next=invector(frontedge2(:,1),...
                 frontedge1(:,[1 2]),ccw(1));
   while last==next
    v1b=v1; v1=mod(v1,numverts1)+1;
    next=invector(frontedge2(:,1),...
                  frontedge1(:,[v1b v1]),ccw(1));
   end
   first=v1;
  end
 
  % a=[9 -9; -3 0 ] , b=[2 -8 4 ; 9 4 -8]
  % was the input that gave the bug...
 
  % Find the second such alternation...
  last=next;
  v1b=v1; v1=mod(v1,numverts1)+1;
  next=invector(frontedge2(:,1),frontedge1(:,[v1b v1]),ccw(1));
  while last==next
   v1b=v1; v1=mod(v1,numverts1)+1;
   next=invector(frontedge2(:,1),frontedge1(:,[v1b v1]),ccw(1));
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
 
  area=abs(det([frontedge2(:,1),...
                supp1(:,hull1(finish))-supp1(:,hull1(start))]));
  return
 end
else % The hardest part is when both Newton polygons 
     % have positive area... 
 % Get ready to search for an initial vert(2)
 v2=1; v2b=numverts2; 

 % Let's find the orientations of the Newton polygons... 
 ccw=[sign(det(frontedge1(:,1:2))) , ... 
      sign(det(frontedge2(:,1:2)))] ; 
 move2=ccw(1)*ccw(2);

 % Now let's find a v2 on hull2 whose angle 
 % cone is disjoint from the angle cone of 
 % the vertex 1 on hull1...  
 % (The main test is to see if a frontedgei(vert(i)+1)  
 % (or -frontedgei(vert(i)) of polygon i is within the 
 % cone generated by frontedgej(vert(j)+1) and 
 % -frontedgej(vert(j)).) There are four ``AND'' 
 % conditions to be ``OR'ed''...
 % NOTE: Sign of polygon j determines the factor in the inequality...
 while conesect(frontedge1(:,[numverts1 1]),ccw(1),...
                frontedge2(:,[v2b v2]),ccw(2)) 
    v2=mod(v2,numverts2)+1 ;  
    v2b=mod(v2b,numverts2)+1 ; 
 end  
end

opt1=1; opt1b=numverts1; 
opt2=v2; opt2b=mod(opt2-2,numverts2)+1;

area=0; % Get ready to compute the mixed area... 
parallelo=zeros(2,4); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's finally define the convex hull of the Minkowski
%% sum...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
move=zeros(2,2); % move(``half'' , polygon ) 
if move2==1             % Depending on the orientations,
 edgeind2=[opt2b opt2]; % the indices of the starting edges of 
 move=[1 -1 ; -1 1];    % the ``forward'' and ``backward'' 
else                    % parallelograms will be different...
 edgeind2=[opt2 opt2b]; % (We'll also need to know which way the 
 move=[1 1 ; -1 -1];    % indices move...)
end
edgeind1=[opt1 opt1b]; 

% Starting from the pinch, incrementally draw the 
% parallelograms, first on one side, and then on the other 
% side, of the first Newton polygon.
for i=1:2
 % Let's first find the forward (or backward) 
 % limit of the edges on the second 
 % Newton polygon, from the pinch...
 lim=0; 
 e1=edgeind1(i); 
 initedge1=ccw(1)*move(i,1)*frontedge1(:,e1);
 e2=edgeind2(i); 
 while det([move(i,2)*frontedge2(:,e2),...
            move(i,1)*initedge1])>0;
  e2=mod(e2+move(i,2)-1,numverts2)+1;
  lim=lim+1; 
 end 

 % Now that we have our limits, let's start 
 % incrementally drawing our parallelograms...
 v1=opt1; edge1=initedge1;
 initedge2=move(i,2)*frontedge2(:,edgeind2(i));
 finaledge2=move(i,2)*frontedge2(:,mod(edgeind2(i)+...
                                   lim*move(i,2)-1,...
                                   numverts2)+1); 
 while det([initedge2,move(i,1)*edge1])>0 
 
  % Do we need to lower the limit?...
  while det([finaledge2,move(i,1)*edge1])<=0
   lim=lim-1; 
   finaledge2=move(i,2)*frontedge2(:,mod(edgeind2(i)+...
                                     lim*move(i,2)-1,...
                                     numverts2)+1); 
  end 
  if lim>=0 
   e2=edgeind2(i); edge2=initedge2; 
   parallelo(1,:)=supp1(1,hull1(v1));
   parallelo(2,:)=supp1(2,hull1(v1));
   for j=0:lim
    area=area+det([edge2,move(i,1)*edge1]);  
    parallelo(:,2)=parallelo(:,1)+edge1; 
    parallelo(:,3:4)=[parallelo(:,2),parallelo(:,4)]+...
                     [edge2,edge2]; 
    parallelo(:,1)=parallelo(:,1)+edge2; 
    e2=mod(e2+move(i,2)-1,numverts2)+1; 
    edge2=move(i,2)*frontedge2(:,e2); 
   end
  end
  v1=mod(v1+move(i,1)-1,numverts1)+1;  
  e1=mod(e1+move(i,1)-1,numverts1)+1;
  edge1=ccw(1)*move(i,1)*frontedge1(:,e1);
 end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here are some functions I use in marea...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This one's for error-handling in the input...
function newsupp=check(supp,ind);  
s=size(supp,2) ; 

if size(supp,1)~=2 
 disp(sprintf(['\nYour ', card(ind),...
               ' support has the wrong size!\n\n']))
 newsupp=[];
 return
end

small=0; 
if s==1 
 disp(sprintf(['\nYour ', card(ind),... 
               ' support has only one point!\n\n']))
 small=1; 
else 
 % (It seems that in our (LL, 6/25/97) applications, 
 %  no support will ever be bigger than, say, 50. 
 %  So we can use the simplest O(n^2) algorithm safely.) 
 repeats=0; i=1; 
 while i<s 
  j=i+1; 
  while j<=s 
   if supp(:,j)==supp(:,i) 
    if repeats==0
     disp(sprintf(['\nYou have repeats in your ', card(ind),...
                   ' support!']))  
     disp(sprintf('(But I''ll remove them because I''m nice...)\n\n'))
     repeats=1; 
    end
    supp=supp(:,[1:j-1 j+1:s]); 
    j=j-1; 
    s=s-1; 
   end
   j=j+1; 
  end 
  i=i+1; 
 end
end
% Check again, just in case.
s=size(supp,2) ; 
if and(s==1,not(small)) 
 disp(sprintf(['\nYour ', card(ind),...
               ' support has only one point!\n\n']))
end
newsupp=supp; 

% This one's for cardinalities in English...
function word=card(ind);
if ind==1
 word='first';
elseif ind==2
 word='second';
elseif ind==3
 word='third'
else 
 word='nth?'
end

% This one's for setting up the Newton polygon 
% (including the edge list) of a support... 
function [numverts,frontedge,hull]=polygon(supp) 
s=size(supp,2); 
if s<=1 
 numverts=s; hull=[]; frontedge=[]; 
 return 
elseif s==2
 hull=[1 2]; numverts=2; 
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
  hull=[amin amax]; 
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
  true=[] ; 
  for i=1:pumverts 
   if det([grontedge(:,mod(i-2,pumverts)+1), ...
           grontedge(:,mod(i-1,pumverts)+1)])~=0 
    true=[true mod(i-1,pumverts)+1 ] ; 
   end
  end 
  % ...and modify as necessary.
  numverts=size(true,2);  
  if numverts<pumverts 
   % Throw out points on edge interiors...
   hull=zeros(1,numverts); 
   hull=gull(true(:)); 
   % Compress the edgelist...
   frontedge=zeros(2,numverts); 
   frontedge(1,:)=supp(1,hull([2:numverts 1]))-...
                   supp(1,hull(1:numverts));  
   frontedge(2,:)=supp(2,hull([2:numverts 1]))-...
                   supp(2,hull(1:numverts));  
  else 
   frontedge=grontedge; 
   hull=gull(1:numverts);
  end
 end
end 

% This is to see if an edge vector lies in the  
% angle cone defined by two consecutive edges...
% (Is inedge of polygon i within the cone generated by 
% -bedge and fedge?)  
% NOTE: The orientation of the polygon that bedge and fedge come 
%       from determines the sign factor in the 
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
 

