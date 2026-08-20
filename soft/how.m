function p=how(a); 

figure(1); 
subplot(1,1,1);
hold off; 
newplot; hold on; 
axis equal; 
axis([-5 5 -5 5 ]) ; 

plot(0,0,'kd'); 

h=a*[0 0 0 1 1 1; 0 1 .5 .5 1 0];  
o=a*[2 2 3 3 2; 0 1 1 0 0];  
w=a*[4 4 4.5 4.5 4.5 5 5; 1 0 0 1 0 0 1 ];  
d=a*[6 6 7 6; 0 1 .5 0];  
y=a*[8 8.5 9 8.5 8.5 ; 1 .5 1 .5 0 ];  

plot(h(1,:),h(2,:));  
plot(o(1,:),o(2,:));  
plot(w(1,:),w(2,:));  
plot(d(1,:),d(2,:));  
plot(y(1,:),y(2,:));  


