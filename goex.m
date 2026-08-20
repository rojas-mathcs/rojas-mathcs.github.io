% this is the tableau i stated in lecture on Tuesday, Jan. 24, 2017
lp1=[20  16  12   0  0  0 -10;
     1   0   0   1  0  0   4;
     2   1   1   0  1  0  10;
     2   2   1   0  0  1  16];  

% the rows of the following matrix are the indices of the 
% bases underlying each step we did for the simplex method. 
bases=[4 5 6; % 10 
       1 5 6; % 90 
       1 2 6; % 122 
       1 2 4; % 146
       2 3 4]; % 154 
% the extra numbers to the right are the value of the 
% objective function at the corresponding bfs... 


% now comes the 2nd tableau from Tuesday, Jan. 24       
lp2=[1 0 0 0 ; [zeros(3,1) lp1([2:4],bases(2,:))]]\lp1;
lp2(1,:)=lp2(1,:)-lp2(1,bases(2,:))*lp2([2:4],:); 

% the 3rd tableau from Tuesday, Jan. 24       
lp3=[1 0 0 0 ; [zeros(3,1) lp2([2:4],bases(3,:))]]\lp2;
lp3(1,:)=lp3(1,:)-lp3(1,bases(3,:))*lp3([2:4],:); 

% the 4th tableau from Tuesday, Jan. 24       
lp4=[1 0 0 0 ; [zeros(3,1) lp3([2:4],bases(4,:))]]\lp3;
lp4(1,:)=lp4(1,:)-lp4(1,bases(4,:))*lp4([2:4],:); 

% the 5th tableau from Tuesday, Jan. 24       
lp5=[1 0 0 0 ; [zeros(3,1) lp4([2:4],bases(5,:))]]\lp4;
lp5(1,:)=lp5(1,:)-lp5(1,bases(5,:))*lp5([2:4],:); 

% these are tableau that i introduced on Thursday, Jan. 26
lp22=lp2; lp22(4,7)=4; 
lp222=lp2; lp222(2,4)=-1; lp222(3,4)=0;  
lp33=[1 0 0 0 ; [zeros(3,1) lp22([2:4],bases(3,:))]]\lp22;
lp33(1,:)=lp33(1,:)-lp33(1,bases(3,:))*lp33([2:4],:); 
lp44=[1 0 0 0 ; [zeros(3,1) lp33([2:4],bases(4,:))]]\lp33;
lp44(1,:)=lp44(1,:)-lp44(1,bases(4,:))*lp44([2:4],:); 
% in reality, i only showed lp22, lp33, lp44, in order to 
% illustrate how a degenerate bfs can lead to the objective 
% function not growing. i was going to also speak of unboundedness, 
% using tableau lp222, but didn't get to. i'll likely speak about 
% lp222 next Tuesday, Jan. 31, 2017... 
