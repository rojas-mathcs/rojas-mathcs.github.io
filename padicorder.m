function padicorder = padicorder(p, a)

if (p < 2)
    disp('P must be greater than 1!');
    padicorder = [];
    return;
end

l = length(a);
padicorder = zeros(size(a));
for i = 1 : l
    padic = 0;
    t = 1;
    k = 0;

    aa = abs(a(i));
    while t <= aa
        t = t*p;
        k = k + 1;
        if (mod(aa, t) == 0)
            padic = k;
        end
    end
    
    padicorder(i) = padic;
end