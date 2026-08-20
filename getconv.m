function [status,conv] = getconv(filename)

% Getconv get data from an OFF file in OOGL format, and returns the data in
% convex hull structure. It can be used as:
% 
%     [status, conv] = getconv(filename)
%
% filename is a string bringing the name of file which contains all the
% data, it is required to be an OOGL file in OFF format. 
% conv returns the data of the geometry object in convex polytope structure
% and if the data is gotten succesfully, status returns as 0;

status = 0;

convfile = fopen(filename, 'r');

tline = fgetl(convfile);
if (lower(tline) ~= 'off')
    disp(strcat(filename,' is  NOT an OFF file !!'));
    status = 2;
    return;
end

% Nv and Nf are the numbers of vertexes and edges in polygon1
temp = fscanf(convfile,'%d %d %d',3);
conv.Nv = temp(1);
conv.Nf = temp(2);
conv.Ne = 0;
conv.vert = [];

%Vert are the coordinate of vertexes
%k = conv.Nv;
%for i = 1 : k
%    tline = fgetl(convfile);
%    while (length(tline) == 0)
%        tline = fgetl(convfile);
%    end
%    temp = sscanf(tline, '%e ');
%    conv.vert = [conv.vert; temp'];
%    keyboard;
%end
    
[conv.vert,count1] = fscanf(convfile, '%e %e %e', [3 conv(1).Nv]);
conv.vert = (conv.vert)';


%edges and facets
k = conv.Nf;
conv.ed = [];
for i = 1: k
    
    tline = fgetl(convfile);
    temp = str2num(tline);
    while (length(temp) == 0)
        tline = fgetl(convfile);
        temp = str2num(tline);
    end
   
    %Nvert = fscanf(convfile, '%d ', 1);
    %temp = fscanf(convfile, '%e ', [1 Nvert+1]);
    %%%
    
    Nvert = temp(1);
    temp = temp(2:Nvert+1);
    
    temp = temp + 1;
    
    conv.fac(i).Nv = Nvert;
    conv.fac(i).v = temp(1:Nvert);
    for j = 1:Nvert
        
        if (j == Nvert)
            Newed1 = [temp(j),temp(1)];
            Newed2 = [temp(1),temp(j)];
        else
            Newed1 = [temp(j), temp(j+1)];
            Newed2 = [temp(j+1), temp(j)];
        end
        
        isolded = 0;
        for m = 1:conv.Ne
            olded = conv.ed(m,:);
            if (compare(olded,Newed1) || compare(olded, Newed2))
                isolded = 1;
                break;
            end
        end
        if (~ isolded) || (conv.Ne == 0)
            conv.ed = [conv.ed; Newed1];
            conv.Ne = conv.Ne + 1;
        end        
    end
end

fclose(convfile);
