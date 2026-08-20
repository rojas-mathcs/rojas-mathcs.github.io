% june 24, 2005 version: 
% todo: fix any bugs which may pop up. 
% questions and comments welcome! 
% J. Maurice Rojas; June 24, 2005
% 
% usage: reamoeba(const,linpow,piclim,parlim,res,fig,flipin,flea) 

function reamoeba(const,linpow,piclim,parlim,res,fig,flipin,flea) 

% const  = cell array of 2xm matrices for multiples:  a1^b1 a2^b2 ... 
% linpow = cell array of 3xm matrices for (a_i t + b_i)^c_i 
% piclim = [minx,maxx,miny,maxy], for picture
% parlim = [minx,maxx] for parameter line 
% res    = density of points on parameter line 
% fig    = which figure for the picture  
% flipin = indices of which signs matter (should be automated later!) 
% flea   = extra sign flips from DFS formula, as arising from Rusek's 
%          implementation... (should be automated later!)  

% # of steps for real part is proportional 
%  to largest of {width,height}...

figure(fig);
hold off; newplot; axis equal; axis xy; axis off; hold on; 

% set the cross...
w=piclim(2)-piclim(1); h=piclim(4)-piclim(3); 
plot([-w w],[0 0],'k-'); plot([0 0],[-h h],'k-');

% pick out which coordinates matter for signs...
% this should be automated later!!! 

tot=(1+parlim(2)-parlim(1))*res 

% get constants and sizes for linlog expressions...
for i=1:2 
 cx(i)=const{i}(2,:)*(log(const{i}(1,:)))';  
end; 

for t=parlim(1):(parlim(2)-parlim(1))/res:parlim(2) 

 for i=1:2
  in{i}=t*linpow{i}(1,:)+linpow{i}(2,:);  
  x(i)=cx(i)+linpow{i}(3,:)*(log(abs(in{i})))'; 
 end; 

 if (x(1)>=piclim(1))&(x(1)<=piclim(2))&(x(2)>=piclim(3))&(x(2)<=piclim(4))
  % note: for some d, one should insert minus signs below!  
  % e.g., 
  % d=3  =>  -+
  % d=54 =>  --
  flipx=flea(1)*prod(sign(in{1}(flipin{1}))); 
  flipy=flea(2)*prod(sign(in{2}(flipin{2}))); 
  % plot(x(1),x(2),'.','Color',[ct/tot 0 1-(ct/tot)]);
  plot(flipx*(x(1)-piclim(1)),flipy*(x(2)-piclim(3)),'b.'); 
  plot(flipy*(x(2)-piclim(3)),flipx*(x(1)-piclim(1)),'b.'); 
 end; 
end; 

 
axis equal;  
title('Amoebified Discriminant Locus for {[6 0 0; 0 3 1],[0 3 1; 6 0 0]}');
xlabel(sprintf('Truncation at [%.1f,%.1f]x[%.1f,%.1f], signs respected...\n Plot obtained via DFS parameteric formula...',piclim(1),piclim(2),piclim(3),piclim(4)));
axis on;

