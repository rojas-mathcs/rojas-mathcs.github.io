function sub2file(convh,numsub,subdiv,alpha)

% sub2file
% sub2file creates an OOGL file from information of a mixed subdivision
% of some polyhodrons.You will be asked to give a file name to save it. Now
% it is valid for two, three or four polyhodrons.sub2file.m uses sub2file2p,
% sub2file3p or sub2file4p to do the job.
%
%   sub2file(convh, numsub, subdiv, alpha)
%   sub2file(convh, numsub, subdiv)
%
% convh is a convex polyhodron stucture array which contains information 
% of the polyhodrons.
% numsub is the number of cells of the subdivion, and subdiv contain the
% infomation for every cell, it give us the whole figure of the subivision
% when combined with convh.
% subdiv is a matrix with each row representing a cell. Every row is orginazed as:
%    [type inex]
% also see sub2file2p, sub2file3p and sub2file4p for more detail.
% alpha controls the ration to shink the cell, it rangs from 0.01 to 1,
% when omitted from the input, the default value 1.0 will be used.
%
% Also see: sub2file2p, sub2file3p and sub2file4p

if nargin < 4
    alpha = 1.0;
elseif nargin < 3
    disp('Not sufficient arguments!!');
    return;
end

if (alpha < 0.01 | alpha > 1.0)
    disp('alpha must be between 0.01 to 1.0 !');
    return;
end

numpoly = length(convh);

if (ndims(subdiv) ~= 2) | (numsub ~= size(subdiv,1)) | (size(subdiv,2) ~= 2*numpoly)
    disp(' Bad subdivision matrix or number !!');
    return;
end

switch numpoly
    case 2
        sub2file2p(convh,numsub,subdiv,alpha);
    case 3
        sub2file3p(convh,numsub,subdiv,alpha);
    case 4
        sub2file4p(convh,numsub,subdiv,alpha);
end

