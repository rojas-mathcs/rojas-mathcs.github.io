function [newa, newc] = polytran2v(a,c)

% polytran2v computes a new polynomial of two variables by replacing the
% original x, y with (x+1) and (y+1). It is used as:
%
%      [newa, newc] = polytran2v(a,c)
%
% a, c is the exponants and coefficients in the old polynomial respectively.
% They have same number of rows, but a has 2 columns and c has 1. 
% newa, newc is the exponants and coefficients of the new polynomail. They
% have not been sorted, but zero items will be kicked out.

m = length(a);

newa = [];
newc = [];
nitems = 0;

for k = 1 : m
    if (c(k) ~= 0)
        
        xorder = a(k,1);
        yorder = a(k,2);
        
        px = ones(1,xorder+1);
        py = ones(1,yorder+1);
        
        for i = 2 : xorder
            px(i) = nchoosek(xorder,i-1);
        end
        
        for i = 2 : yorder
            py(i) = nchoosek(yorder,i-1);
        end
        
        for i = 1 : xorder + 1
            for j = 1 : yorder + 1
                temporder = [i-1,j-1];
                tempc = px(i) * py(j) * c(k);
                hasadded = 0;
                for m = 1 : nitems
                    if compare(temporder, newa(m,:))
                        newc(m) = newc(m) + tempc;
                        hasadded = 1;
                        break;
                    end
                end
                
                if ~hasadded
                    newa = [newa; temporder];
                    newc = [newc; tempc];
                    nitems = nitems + 1;
                end
            end
        end
    end
end

k = 1;
while k <= nitems
    if newc(k) == 0
        if nitems == 1
            nitems = 0;
            newc = [];
            newa = [];
        else
            if k == 1
                newa = newa(2:nitems, :);
                newc = newc(2:nitems);
            elseif k == nitems
                newa = newa(1:nitems-1,:);
                newc = newc(1:nitems-1);
            else
                newa = [newa(1:k-1,:); newa(k+1:nitems,:)];
                newc = [newc(1:k-1); newc(k+1:nitems)];
            end
        end
        nitems = nitems - 1;
    else
        k = k + 1;
    end
end
        