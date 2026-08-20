function [lfacbak, convh, eepic] = pamoeba(a, c, p)

% pamoeba computes the amoeba and show it in figure. Now, it is valid for
% two variables polynomial. a is a 2 columns matrix reprsenting the
% exponant of items in the polynomial and c is the corresponding coefficients.
% p is a prime value for p-adic order.
% the function can be called as:
%
%        [lfacbak, convh, eepic] = pamoeba(a, c, p)
%
% for the returned value: 
% lfacbak is an array of stucture containg information about the lower
% facets of the convex hull
% convh is a structure type variable containing information about the convex
% hull
% eepic is a stucture type variable discribing the figure, including the
% points, arrows and lines.
%
x = [a, padicorder(3,c)];

    [convh, kface] = convhullng(x);
    
    w = [];
    centconv = zeros(1,3);
    
    for j = 1 : convh.Nv
        centconv = centconv + convh.vert(j,:);
    end
    centconv = centconv / convh.Nv;
    
    lfacbak = [];
    nlowf = 0;
    
    eepic.Nv = 0;
    eepic.Na = 0;
    eepic.Nl = 0;
    eepic.v = [];
    eepic.a = [];
    eepic.l = [];
    
    hold on;
    axis equal;
    
    for fi = 1 : convh.Nf
        
        centface = centeval(2,1,fi,convh,centconv);
        facvec = [];
        for j = 1 : convh.fac(fi).Nv
            p1 = convh.fac(fi).v(j);
            facvec = [facvec; convh.vert(p1,:) - centface];
        end
        
        facenorm = null(facvec);
        if abs(facenorm(3)) < 1.0e-10 
            islowerf = 0;
        else
            facenorm = facenorm / facenorm(3);
            if (centconv - centface) * facenorm > 0
                islowerf = 1;
            else
                islowerf = 0;
            end
        end
            
        if islowerf
            
            %draw the point
            newlowf.fi = fi;
            newlowf.norm = facenorm';
            newlowf.cent = centface;
            
            lfacbak = [lfacbak; newlowf];
            nlowf = nlowf + 1;
            
            plot(facenorm(1), facenorm(2),'k.');
            eepic.Nv = eepic.Nv + 1;
            eepic.v = [eepic.v; [facenorm(1), facenorm(2)]];
            
            for ei = 1 : convh.fac(fi).Nv
                if ei == convh.fac(fi).Nv
                    p1 = convh.fac(fi).Nv;
                    p2 = 1;
                else
                    p1 = ei;
                    p2 = ei + 1;
                end
                p1 = convh.fac(fi).v(p1);
                p2 = convh.fac(fi).v(p2);
                
                tested = [p1, p2];
                edvec = convh.vert(p1,:) - convh.vert(p2,:);
                cented = (convh.vert(p1, :) + convh.vert(p2,:)) / 2;
                edvec = edvec(1:2);
                cented = cented(1:2);
                ednorm = null(edvec);
                if (centface(1:2) - cented) * ednorm < 0
                    ednorm = -1 * ednorm;
                end
                
                ff = (a * ednorm)';
                aa = eye(convh.Nv) * -1.;
                bb = zeros(convh.Nv,1);
                aeq = ones(1,convh.Nv);
                beq = 1;
                lb = zeros(convh.Nv,1);
                ub = ones(convh.Nv,1);
                
                [lamda, h] = linprog(ff,aa,bb,aeq,beq,lb,ub);
                
                if (h - cented*ednorm) > -1.0e-9
                    xyline = [facenorm(1:2)'; facenorm(1:2)' + ednorm'];
                    plot(xyline(:,1),xyline(:,2),'k');
                    eepic.Na = eepic.Na + 1;
                    eepic.a = [eepic.a; [facenorm(1:2)', ednorm']];
                end
            end
            
            for lfi = 1 : nlowf - 1
                %Is they share the same edge?
                commonedge = 0;
                
                for i = 1 : convh.fac(fi).Nv
                    if i == convh.fac(fi).Nv
                        p1 = convh.fac(fi).Nv;
                        p2 = 1;
                    else
                        p1 = i;
                        p2 = i + 1;
                    end
                    newed = [convh.fac(fi).v(p1), convh.fac(fi).v(p2)];
                    
                    oldfi = lfacbak(lfi).fi;
                    for j = 1 : convh.fac(oldfi).Nv
                        if j == convh.fac(oldfi).Nv
                            p1 = convh.fac(oldfi).Nv;
                            p2 = 1;
                        else
                            p1 = j;
                            p2 = j + 1;
                        end
                        olded = [convh.fac(oldfi).v(p1), convh.fac(oldfi).v(p2)];
                        olded1 = [convh.fac(oldfi).v(p2), convh.fac(oldfi).v(p1)];
                        if compare(olded, newed) | compare(olded1, newed)
                            commonedge = 1;
                            break;
                        end
                    end
                    
                    if commonedge 
                        break;
                    end
                end
                        
                if commonedge
                    xyline = [facenorm(1:2)'; lfacbak(lfi).norm(1:2)];
                    plot(xyline(:,1),xyline(:,2),'k');
                    eepic.Nl = eepic.Nl + 1;
                    eepic.l = [eepic.l; [facenorm(1:2)', lfacbak(lfi).norm(1:2)]];
                end
            end
        end
    end
     