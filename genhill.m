% copyright 2012/09/09, j. maurice rojas 
% This code encrypts plain-text as I requested in HW#1 of 
% Math 470 (Fall 2012). 
% 
% usage: cipher=genhill(plain,a,b) 
% where...
%  a = 2x2 integer matrix with determinant relatively prime to 26 
%  b = a 1x2 vector of integers 
%  plain = a string 
%  cipher = the string obtained by replacing plain by a 
%           a vector V with elements in Z/26Z (in the usual way with 
%           a->0, b->1, ..., z->25), hitting consecutive pairs of elements 
%           of V with the map x -> b+xa, then converting back into text. 
%
% note1: the plain-text is assumed to consist solely of lower-case letters,   
%        with NO punctuation, numbers, or spaces. 
% note2: if the plain text is of odd length then the last character is 
%        IGNORED.  
function cipher=genhill(plain,a,b); 

% format rat; 

len=size(plain,2);

cipher=['']; 
for i=1:2:len
 x=double(plain([i i+1]))-97;   
 y=mod(b+x*a,26); 
 cipher=[cipher char(y+97)];
end; 

len 


