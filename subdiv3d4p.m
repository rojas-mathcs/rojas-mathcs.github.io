function [k,subdiv] = subdiv3d4p(conv0,a)

%
% subdiv3d4p
%Subdiv3d4p is a function to find the mixed subdivision of four convex polyhodrons in 
%3D space. conv0 is a convex polyhoron array which brings the information
%of the four polyhodrons, and a is a 3 by 4 matrix used for lifting.
%          
%         [k, subdiv] = subdiv3d4p(conv0, a)
%
%         [k, subdiv] = subdiv3d4p(conv0)
%
%For the return value:
%   k is the number of cells
%   subdiv is k by 8 matrix, every row for a cell. And every row is divided
%   to two parts:
%         [type | index]
%   The type can be [0 0 0 3], [0 0 1 2], [ 1 0 1 1] etc.
%
% If a is omitted in the input, as:
%         [k, subdiv] = sbudiv3d4p(conv0)
% it is creted automatically.


outfile = fopen('outfile', 'w');
convg = conv0;

% put them in 3d space and linearly lift polygon2
if (nargin < 2)
    a = 1000. * rand(3, 4);
    a(:,1) = zeros(3, 1);
end

for i = 1 : 4
    convg(i).vert = [convg(i).vert, convg(i).vert * a(:,i)];
end

% these variable are for linear programing
V = [(convg(1).vert)',(convg(2).vert)',(convg(3).vert)',(convg(4).vert)'];
f = [0 0 0 1] * V;

totalNv = 0;
for i = 1:4
    totalNv = totalNv + convg(i).Nv;
end

Aeq = [ones(1, convg(1).Nv), zeros(1, convg(2).Nv+convg(3).Nv+convg(4).Nv)];
Aeq = [Aeq; zeros(1, convg(1).Nv), ones(1, convg(2).Nv), zeros(1,convg(3).Nv+convg(4).Nv)];
Aeq = [Aeq; zeros(1, convg(1).Nv + convg(2).Nv), ones(1, convg(3).Nv),zeros(1,convg(4).Nv)];
Aeq = [Aeq; zeros(1, convg(1).Nv + convg(2).Nv + convg(3).Nv), ones(1, convg(4).Nv)];
Aeq = [Aeq; V(1:3,:)];

A = eye(totalNv) * -1.;
b = zeros(totalNv, 1);
lb = zeros(totalNv, 1);
ub = ones(totalNv, 1);

% centv1 and centv2 are the central point in polygon1 and polygon2, they
% are used in type (0 2) and (2 0) to reduce computation burdun.
centconv = zeros(4,4);
for i = 1: 4
    for j = 1 : convg(i).Nv
        centconv(i, :) = centconv(i, :) + convg(i).vert(j,:);
        
    end
    centconv(i, :) = centconv(i, :) / convg(i).Nv;
end

k = 0;
subdiv = [];

for t1 = 0:3
    for t2 = 0 : 3
        for t3 = 0:3
            for t4 = 0:3
                if (t1+t2+t3+t4 == 3) 
                    k1 = getk(t1,1,convg);
                    k2 = getk(t2,2,convg);
                    k3 = getk(t3,3,convg);
                    k4 = getk(t4,4,convg);
                    for i1 = 1 : k1
                        for i2 = 1 : k2
                            for i3 = 1 : k3
                                for i4 = 1 : k4
                                    cent0 = 0;
                                    cent0 = cent0 + centeval(t1,1,i1,convg,centconv);
                                    cent0 = cent0 + centeval(t2,2,i2,convg,centconv);
                                    cent0 = cent0 + centeval(t3,3,i3,convg,centconv);
                                    cent0 = cent0 + centeval(t4,4,i4,convg,centconv);
                                    beq = [1; 1; 1;1;(cent0(1:3))'];
                                    [lamda, h] = linprog(f,A,b,Aeq,beq,lb,ub);
                                    if (h >= cent0(4))
                                        k = k + 1;
                                        [t1 t2 t3 t4 k1 k2 k3 k4 i1 i2 i3 i4]
                                        lamda
                                        h
                                        cent0
                                        subdiv = [subdiv; t1 t2 t3 t4 i1 i2 i3 i4];
                                        fprintf(outfile, '%d %d %d %d %d %d\n', t1, t2, t3, t4, i1, i2, i3, i4);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

fclose(outfile);
