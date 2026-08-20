function [k,subdiv] = subdiv3d(conv0, a)

%
% subdiv3d
% Subdiv3d2p is a function to find the mixed subdivision of some convex polyhodrons in 
%3D space. Now it is valid for two, three or four convex polyhodrons. 
%
%conv0 is a convex polyhoron array which brings the information
%of the polyhodrons, and a is a matrix used for lifting.
%          
%         [k, subdiv] = subdiv3d2p(conv0, a)
%
%For the return value:
%   k is the number of cells
%   subdiv is k a matrix. Every row reprents a cell and is divided into two
%   parts:
%         [type | index]
%
%subdiv3d uses subdiv3d2p, subdiv3d3p or subdiv3d4p to do the job. 
%
% If a is omitted in the input, as:
%         [k, subdiv] = sbudiv3d2p(conv0)
% it is creted automatically.
%
% Also see: subdiv3d2p, subdiv3d3p and subdiv3d4p


numpoly = length(conv0);

if (nargin == 2)
    if (ndims(a) ~= 2)
        disp('The lifting matrix must be 2D !!')
        return;
    end
    
    [m, n] = size(a);
    if n ~= numpoly
        disp('The number of columns of the lifting matrix must correspond with the number of polyhodrons !!');
        return;
    elseif m ~= 3
        disp('The lifting matrix need three rows!!');
        return;
    end
    
end

if (nargin < 2)
    switch numpoly
        case 2
            [k, subdiv] = subdiv3d2p(conv0);
        case 3
            [k, subdiv] = subdiv3d3p(conv0);
        case 4
            [k, subdiv] = subdiv3d4p(conv0);
    end
else
    switch numpoly
        case 2
            [k, subdiv] = subdiv3d2p(conv0, a);
        case 3
            [k, subdiv] = subdiv3d3p(conv0, a);
        case 4
            [k, subdiv] = subdiv3d4p(conv0, a);
    end
end


        
    
    