% grades for whole course...
hha=[ 
% 5  5  0 20  0  0; is  michael sid alexander dead? 
10 10  0 20 63 10;  
 5 10 10 10 44  0;  
10  5 10 20 42  0; 
10 10 15 15 87 10;  
 5 10 10 10 55  0;  
10  0 15  0 67  0;  
10 10 10 20 80 10; 
10  5 15 20 47  0;  
 5 10  5 10 25  2;  
10 10  5 15 20  2; 
 5 10 15 11 49  0;  
10 10  5 15 63  0;  
% 0  0  0  0  0  0; is derek jason weiner dead?
 5  5 10  0 40  2];  

% (key: q1/10 hw/10 q2/10 q3/20 mt/87 ec_for_mt/10) 
ha=hha*[.8;.8;.8;.8;6/9.7;6/9.7];
means=mean(ha); 

sds=std(ha); 

mins=min(ha); 
maxs=max(ha); 
meds=median(ha); 

hist(ha,10) ; 
title('MATH 302 (10:00-11:25): Score Distribution for First Half (Summer 2003)');
xlabel(sprintf('MAX=%d , HI=%d  ,  MEAN+-SD = %.1f+-%.1f  ,  MEDIAN=%.1f  ,  LO=%.1f\n%d Students',100,maxs,means,sds,meds,mins,size(ha,1))); 

