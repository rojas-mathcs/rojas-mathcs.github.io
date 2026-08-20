% 3/17/09, copyright J. Maurice Rojas
% This function sets up the string for the plot title for 
% the drawings done by chambers.m . 
% Usage: 
% polystring=setdiscplottitle(a,n,m,cell,cano,dots); 
% where...
%  a       = support of the family of polynomials you're working with
%  m       = # of support points 
%  n       = # of variables in your family of polynomials  
%  cell    = n-subset of [1,...,m] with det(a(:,cell)) odd...
%  cano    = 1 or 0, according as you're about to plot a canonical slice or not
%  dots    = 1 or 0 according as you want some random polynomials drawn or not
% and the output is the string that will be displayed as the title of 
% your plot...
function polystring=setdiscplottitle(a,n,m,cell,cano,dots)  

origif=(abs(a(:,1))==zeros(n,1));
if dots==1
 polystring='%d real centered Gaussian polynomials, log(variance)=[%1.2f %1.2f %1.2f %1.2f'; 
 for j=5:m
  polystring=[polystring,' %d'];
 end;
 if cano==0
  polystring=[polystring,'],\n\n scaled down to '];
 else
  polystring=[polystring,'],\n\n projected to canonical slice of c_1'];
 end;
 if sum(origif(:))<n
  polystring=[polystring,'x^{']; 
  if n==1 
   polystring=[polystring,sprintf('%d}',a(1))]; 
  else 
   polystring=[polystring,sprintf('[%d',a(1,1))];
   for i=2:n
    polystring=[polystring,sprintf(' %d',a(i,1))];
   end;
   polystring=[polystring,']}'];
  end; 
 else 
  if cano==0 
   polystring=[polystring,'1']; 
  end; 
 end; 
else 
 if cano==0
  polystring='Slice of Nabla_A(R) plotted on log paper, for the family\n ';
 else
  polystring='Canonical slice of Nabla_A(R), plotted on log paper, for the family\n c_1';
 end;
 if sum(origif(:))<n
  polystring=[polystring,'x^{'];
  if n==1
   polystring=[polystring,sprintf('%d}',a(1))];
  else
   polystring=[polystring,sprintf('[%d',a(1,1))];
   for i=2:n
    polystring=[polystring,sprintf(' %d',a(i,1))];
   end;
   polystring=[polystring,']}'];
  end;
 else
  if cano==0
   polystring=[polystring,'1'];
  end;
 end;
end;

% either way, we need to print the shape of the homogenized family
% of polynomials being presented...
for j=2:m
 if j~=cell
  polystring=[polystring,sprintf('+c_{%d}x^{',j)];
 else
  polystring=[polystring,sprintf('+x^{')];
 end;
 if n==1
  polystring=[polystring,sprintf('%d}',a(j))];
 else
  polystring=[polystring,sprintf('[%d',a(1,j))];
  for i=2:n
   polystring=[polystring,sprintf(' %d',a(i,j))];
  end;
  polystring=[polystring,']}'];
 end;
end;

% continue setting up title...
if dots==1
 if n==1
  if cano==0
   polystring=[polystring,sprintf(',\n plotted atop Nabla_{[%d %d %d %d]}(R) on log paper',a(1),a(2),a(3),a(4))];
  else
   polystring=[polystring,sprintf(',\n plotted atop canonical slice of Nabla_{[%d %d %d %d]}(R) on log paper',a(1),a(2),a(3),a(4))];
  end;
 else
  if cano==0
   polystring=[polystring,sprintf(',\n plotted atop Nabla_A(R) on log paper')]; 
  else
   polystring=[polystring,sprintf(',\n plotted atop canonical slice of Nabla_A(R) on log paper')]; 
  end;
 end;
end;
