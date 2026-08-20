y0=-1.5; % ymin
z0=-2; % zmin
z1=1;  % zmax
n=20;  % grid size - 1 

z=z0+(0:n)*(z1-z0)/n;               % range of z's
t=sqrt((sqrt(z.^2-12*y0)-z)/6);     % t to give correct y0
z=z'*ones(1,n+1);                   % array of z's
t=t'*(-1+(0:n)*2/n);                % array of t's
x=-4*t.^3-2*z.*t;
y=x.*t+z.*t.^2+t.^4;
h=surf(x,y,z);
axis('square');rotate(h,[0 90],180);rotate(h,[0 0],90);
axis([-4 4 -2 2 -2 2]); rotate(h,[0 90],30);

