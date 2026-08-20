% grades for whole course...
hhb=[ 5  5 10   20 68 9; 
 5 10 15   20 68 0;  
 0 10  0   20 35 0;  
 5 10  2.5 10 26 0;  
 5 0   0   15 52 0;   
% 0 0   0    0  0 0; is duc ho kwon dead?   
 5 0   0    0 28 0;  
 5 10  5   10 54 0; 
10 10 15   20 68 2;  
 0  0  5    5  0 0;   
10  5 10   10 45 0;     
 5  5 10    0 61 0;  
10 10  5   10 46 0;  
 5  5  5   15 72 0;  
 5  5  0   10 38 0;   
 0 10  0    0 18 0 ];  

% (key: q1/10 hw/10 q2/10 q3/20 mt/87 ec_for_mt/10) 
hb=hhb*[.8;.8;.8;.8;6/9.7;6/9.7];
means=mean(hb); 

sds=std(hb); 

mins=min(hb); 
maxs=max(hb); 
meds=median(hb); 

hist(hb,10) ; 
title('MATH 302 (12:00-13:25): Score Distribution for First Half (Summer 2003)');
xlabel(sprintf('MAX=%d , HI=%.1f  ,  MEAN+-SD = %.1f+-%.1f  ,  MEDIAN=%.1f  ,  LO=%.1f\n%d Students',100,maxs,means,sds,meds,mins,size(hb,1))); 

