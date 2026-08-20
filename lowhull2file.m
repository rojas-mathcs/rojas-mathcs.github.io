function lowhull2file(lfacbak,convhull,c)

% lowhull2file saves the lower hull to a OOGL file. It is used as 
%
%       lowhull2file(lfacbak, convhull, c)
%
% lfacbak comes from pamoeba.m. It gives the index of all lower faces in
% the convex hull represented by convhull. c is the color used. It is an
% vector including 3 elements as RGB color.
% You will be asked to give the file name.

if (nargin < 5)
    filename = input('Give the name of the file to save the low hull:','s');
end

if (length(filename) == 0)
    filename = 'matlab.off';
end

[ooglfid, message] = fopen(filename, 'w+');
    
if ooglfid < 0
        
    s = strcat('Error in openning file :',message, '\n');
    display(s);
    return;
    
else
    
    numpoints = convhull.Nv
    numfaces = length(lfacbak);
    
    fprintf(ooglfid,'%s\n','OFF');
    fprintf(ooglfid,'%d   %d  %d \n',numpoints, numfaces, 0);
        
    for i = 1 : numpoints
        fprintf(ooglfid, '%-13.6e', convhull.vert(i,:));
        fprintf(ooglfid, '%s\n', ' ');
    end
        
    for i = 1 : numfaces
        ffi = lfacbak(i).fi;
        fprintf(ooglfid, '%d   ', convhull.fac(ffi).Nv);
        fprintf(ooglfid, '%6d  ', convhull.fac(ffi).v - 1);
        fprintf(ooglfid, '%8.3f ', c);
        fprintf(ooglfid, '%s\n', ' 1.0');
        %keyboard
    end
        
    fclose(ooglfid);
end
