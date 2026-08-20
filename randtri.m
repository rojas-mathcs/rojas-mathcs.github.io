a=[0 3 7];
gap1=a(2)-a(1)-1;
gap2=a(3)-a(2)-1;
num=100;

hold off; subplot(1,1,1); newplot; hold on; 

phase2=exp(2*pi*sqrt(-1)*rand); 
phase1=exp(2*pi*sqrt(-1)*rand); 
for i=1:num
 plot(roots([1 zeros(1,gap2) (8+rand)*phase2 zeros(1,gap1) phase2]),'b.');
end; 

axis equal; 
