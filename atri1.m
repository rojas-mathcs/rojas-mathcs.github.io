% EXAMPLE 3.5: 
% 
% The example below is a new  
% attempt at visualizing the discriminant of the 
% polynomial system 
% 
% x^4+ay^2+y=0
% y^4+bx^2+x=0 
%
% The underlying principal A-determinant 
% was calculated via maple.  Theory (looking at 
% the mixed volume bounds for resultant degrees) 
% implies that the principal A-determinant should have 
% degree <= 66.  A Sylvester-in-cascade computation 
% then gave a MULTIPLE --- of degree 144 --- of the 
% principal A-determinant.  This multiple factored 
% into polynomials of degree 36 and 108, so it would 
% appear that our desired discriminant is exactly 
% factor of degree 36, which is presented and 
% amoebified below...
% 
% J. Maurice Rojas, June 22, 2005...
% 

supp=[11 16; 4 14; 6 6; 15 15; 14 19 ; 18 18 ; 7 17 ; 11 1 ;...
      16 11 ; 8 13 ; 9 9 ; 12 12 ; 17 7 ; 3 3 ; 15 0 ; 0 15 ; ... 
      14 4 ; 5 10 ; 1 11 ; 7 2 ; 19 14 ; 13 8 ; 10 5 ; 2 7 ; 0 0 ]';
coeffs=[-1544524446302208 -493055220453376 548511851520000 155866644873216 ...
        35664401793024 -8916100448256  1259052145311744 -58574067632640000 ...
       -1544524446302208 -6670333073424384 5813185762672640 8039729209540608 ...
       1259052145311744 14860665000000000 11112006825558016 11112006825558016 ...
       -493055220453376 1740077249691648 -58574067632640000 50356965330000000 ...
       35664401793024 -6670333073424384 1740077249691648 50356965330000000 ...
       437893890380859375];
lim=[-1,1,-1,1]; nn=10; pp=5; sig=1; fig=1;
