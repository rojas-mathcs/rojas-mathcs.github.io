function [convh, kface] = convhullng(x)

%CONVHULLNG creative the OOGL file in OFF format of the convex hull created
%by convhulln
%
% [convh, kface] = convhullng(x) returns a structure of convex hull convh
%              and the index from convhulln
%
% convhullng uses convhull.m in Matolab
% 
% Filename asked while runing can be any format you like, Geomview can recognize
% it by checking the first line in the output file.
%
% See also convhulln, convhull, qhull
% 
% Writen at December 14, 2002

if nargin < 1
    error('Need argument to be inputed.');
else
    
    if isempty(x) kk=[]; return, end;
    
    if size(x, 2) ~= 3
        error('This programe is used for 3D!')
        return;
    end
    
    filename = input('Give the name of the file to save the convex hull:','s');
    
    [ooglfid, message] = fopen(filename, 'w+');
    
    if ooglfid < 0 
        s = strcat('Error in openning file :',message, '\n');
        error(s);
        return;
    else 
        
        numpoints = size(x,1);
        convh.Nv = numpoints;
        convh.vert = x;
        xx = x - ones(size(x,1),1) * x(1, :);
        
        if numpoints == 1
            convh.Ne = 0;
            convh.ed = [];
            convh.Nf = 0;
            convh.fac = [];
        elseif rank(xx) == 1
            [temp i] = max(max(xx) - min(xx));
            [temp i1] = max(xx(:,i));
            [temp i2] = min(xx(:,i));
            convh.Ne = 1;
            convh.ed = [i1 i2];
            convh.Nf = 0;
            convh.fac = [];
        elseif numpoints == 3
            convh.Ne = 3;
            convh.ed = [1 2; 2 3; 3 1];
            convh.Nf = 1;
            convh.fac(1).Nv = 3;
            convh.fac(1).v = [1 2 3];
        else
            
            kk = convhulln(x);
        
            convh.Nf = 0;
            convh.Ne = 0;
            numedges = 0;
            convh.ed = [];
            kface = kk;
            numfaces = size(kk,1);
            
            while numfaces > 0
                %New face
                v = x(kk(1,2),:) - x(kk(1,1),:);
                v = [v; x(kk(1,3),:) - x(kk(1,1),:)];
                newfac = kk(1,:);
                
                if rank(v) == 1
                    if numfaces == 1
                        kk = [];
                    else
                        kk = kk(2:numfaces, :);
                    end
                    numfaces = size(kk,1);
                    
                elseif numfaces == 1
                    %the only face remained
                    convh.Nf = convh.Nf + 1;
                    convh.fac(convh.Nf).Nv = 3;
                    convh.fac(convh.Nf).v = newfac;
                    numfaces = 0;
                    
                else
                    %find all the face in the same plane
                    p = null(v);
                    newfacid = [1];
                    newed = [newfac(1) newfac(2); newfac(2) newfac(3); newfac(3) newfac(1)];
                    newNe = 3;
                    newv = newfac';
                    newvc = [x(newv(1),:); x(newv(2),:);  x(newv(3),:)];
                    
                    for ifac = 2 : numfaces
                        tfac = kk(ifac,:);
                        tv = x(tfac(2),:) - x(newfac(1), :);
                        tv = [tv; x(tfac(3),:) - x(newfac(1),:)];
                        tv = [tv; x(tfac(1),:) - x(newfac(1),:)];
                        tv = abs(tv * p);
                        
                        if tv(1) < 1.0e-12 & tv(2) < 1.0e-12 & tv(3) < 1.0e-12
                            % a face in the same plane
                            newfac = [newfac; tfac];
                            newfacid = [newfacid; ifac];
                            
                            % new verticle
                            for i = 1 : 3
                                if ~isold(tfac(i),newv)
                                    newv = [newv; tfac(i)];
                                    newvc = [newvc; x(tfac(i),:)];
                                end
                            end
                            
                            % new edges
                            for i = 1 : 3
                                if i == 3
                                    newed1 = [tfac(3), tfac(1)];
                                else
                                    newed1 = [tfac(i), tfac(i+1)];
                                end
                                if ~isold(newed1, newed)
                                    newed = [newed; newed1];
                                    newNe = newNe + 1;
                                end
                            end
                        end    
                    end
                    
                    % find the edges of the comman face
                    gooded = [];
                    ngooded = 0;
                    for ied = 1: newNe
                        edv = x(newed(ied,1), :) - x(newed(ied,2), :);
                        edv = edv / norm(edv,2);
                        edv = [edv; p'];
                        edp = null(edv);
                        mided = (x(newed(ied,1),:) + x(newed(ied,2),:))/2;
                        temp = (newvc - ones(size(newvc,1),1)*mided) * edp;
                        np = sum(temp > 1.0e-12);
                        nn = sum(temp < -1.0e-12);
                        if nn * np == 0
                            if ngooded > 0
                                isoneline = 0;
                                for i = 1 : ngooded
                                    linecheck = [0 0 0; edv(1,:)];
                                    linecheck = [linecheck; x(gooded(i,1),:) - x(newed(ied,2),:)];
                                    linecheck = [linecheck; x(gooded(i,2),:) - x(newed(ied,2),:)];
                                    linecheck1 = linecheck(2:3,:);
                                    linecheck2 = linecheck([2,4],:);
                                    if (rank(linecheck1) == 1 & rank(linecheck2) == 1)
                                        isoneline = 1;
                                        break;
                                    end
                                end
                                if isoneline
                                    linev = [newed(ied,2);newed(ied,1); gooded(i,1); gooded(i,2)];
                                    [temp, i0] = max(max(linecheck) - min(linecheck));
                                    [temp, i1] = max(linecheck(:,i0));
                                    [temp, i2] = min(linecheck(:,i0));
                                    gooded(i,:) = [linev(i1), linev(i2)];
                                else
                                    gooded = [gooded;newed(ied,:)];
                                    ngooded = ngooded + 1;
                                end
                            else
                                gooded = [gooded; newed(ied, :)];
                                ngooded = ngooded + 1;
                            end
                        end
                        %keyboard;
                    end
                    %keyboard;
                    
                    % find the verticals of the comman face
                    % define the new face
                    
                    finalNv = 2;
                    finalv = gooded(1,:);
                    finalNe = 1;
                    finaled = gooded(1,:);
                    
                    ngooded = size(gooded, 1);
                    gooded = gooded(2:ngooded, :);
                    ngooded = ngooded - 1;
                    
                    %keyboard;
                    while ngooded > 1
                        lastp = finalv(finalNv);
                        lastlastp = finalv(finalNv-1);
                        lasts = x(lastp, :) - x(lastlastp, :);
                        foundnext = 0;
                        for i = 1 : ngooded
                            if (gooded(i, 1) == lastp)
                                foundnext = 1;
                                nextp = gooded(i,2);
                                break;
                            elseif (gooded(i,2) == lastp)
                                foundnext = 1;
                                nextp = gooded(i,1);
                                break;
                            end
                        end
                        
                        if nextp == finalv(1), break, end;
                        
                        if ~foundnext
                            keyboard;
                            error('error in construct new face!');
                            return;
                        else
                            finalv = [finalv, nextp];
                            finalNv = finalNv + 1;
                            finaled = [finaled; [lastp, nextp]];
                            finalNe = finalNe + 1;
                        end
                        
                        if ngooded == 1
                            gooded = [];
                        elseif i == ngooded
                            gooded = gooded(1:ngooded-1, :);
                        elseif i == 1
                            gooded = gooded(2:ngooded, :);
                        else
                            gooded = [gooded(1:i-1,:); gooded(i+1 : ngooded,:)];
                        end
                        ngooded = ngooded - 1;
                        %keyboard;
                    end
                    
                    lastp = finalv(finalNv);
                    firstp = finalv(1);
                    lastlastp = finalv(finalNv - 1);
                    
                    tempv = x(finalv(2),:) - x(finalv(1), :);
                    tempv = [tempv; x(firstp,:) - x(lastp, :)];
                    tempv = [tempv; x(lastp,:) - x(lastlastp, :)];
                    vrank = rank(tempv);
                    
                    if (vrank == 1)
                        finalv = finalv(2 : finalNv-1);
                        finalNv = finalNv - 2;
                        finaled(1,1) = lastlastp;
                        finaled = finaled(1:finalNe-1,:);
                        finalNe = finalNe - 1;
                    elseif (rank(tempv(1:2, :)) == 1)
                        finaled(1,1) = lastp;
                        finalv = finalv(2:finalNv);
                        finalNv = finalNv - 1;
                    elseif (rank(tempv(2:3, :)) == 1)
                        finaled(finalNe,2) = firstp;
                        finalv = finalv(1:finalNv - 1);
                        finalNv = finalNv - 1;
                    else
                        finalNe = finalNe + 1;
                        finaled = [finaled; [lastp firstp]];
                    end
                    
                    if (convh.Nf == 0)
                        convh.Nf = 1;
                        convh.Ne = finalNe;
                        convh.ed = finaled;
                        convh.fac(1).Nv = finalNv;
                        convh.fac(1).v = finalv;
                    else
                        convh.Nf = convh.Nf + 1;
                        convh.fac(convh.Nf).Nv = finalNv;
                        convh.fac(convh.Nf).v = finalv;
                        %keyboard;
                        for i = 1 : finalNe
                            newed1 = finaled(i, :);
                            if ~isold(newed1, convh.ed)
                                convh.ed = [convh.ed; newed1];
                                convh.Ne = convh.Ne + 1;
                            end
                        end
                    end
                    
                    % clear faces used in kk, change numfaces
                    nnewfac = length(newfacid);
                    for i = 1 : nnewfac
                        kk(newfacid(i), :) = [0 0 0];
                    end
                    
                    i = 1;
                    nkk = size(kk,1);
                    
                    while (i <= nkk)
                        if compare(kk(i,:), [0 0 0])
                            if nkk == 1
                                kk = [];
                            elseif i == 1
                                kk = kk(2:nkk, :);
                            elseif i == nkk
                                kk = kk(1:nkk-1, :);
                            else
                                kk = [kk(1:i-1,:);kk(i+1:nkk,:)];
                            end
                            nkk = nkk -1;
                        else
                            i = i + 1;
                        end
                        %keyboard
                    end
                    numfaces = size(kk,1);
                    %keyboard
                end               
                
            end
        end
            
        fprintf(ooglfid,'%s\n','OFF');
        fprintf(ooglfid,'%d   %d  %d \n',convh.Nv, convh.Nf, 0);
        
        for i = 1 : convh.Nv
            fprintf(ooglfid, '%-13.2f ', convh.vert(i,:));
            fprintf(ooglfid, '%s\n', ' ');
        end
        
        for i = 1 : convh.Nf
            %colorarray = rand(1,3);
            fprintf(ooglfid, '%6d  ', convh.fac(i).Nv);
            fprintf(ooglfid, '%6d  ', convh.fac(i).v - 1);
            %fprintf(ooglfid, '%d', i);
            %fprintf(ooglfid, '%10.4f',colorarray);
            fprintf(ooglfid, '%s\n', ' ');
        end
        
        fclose(ooglfid);
    end
end
