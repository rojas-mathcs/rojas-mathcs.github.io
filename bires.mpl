with(linalg);
with(geometry);


# supp = input support
# vars = variables
# return the support as a two dimensional matrix


supp_mat := proc (supp::list, vars::list)
  local m, v, nmons, nvars, suppmat;
  

  nmons := nops(supp);
  nvars := nops(vars);
  nvars := nops (vars);
  suppmat := matrix(nmons, nvars);

  for m from 1 to nmons do
    for v from 1 to nvars do suppmat[m,v] := degree (supp[m], vars[v]); od;
  od;	# for monomial

  evalm (suppmat);
end:	# supp_mat

# input: support matrix
# output: list of vertices in counterclockwise order


convhull := proc(suppmat, nmons)

local pts, m, hull, vertlist, v;

for m from 1 to nmons do
pts[m] := point(p||m, convert(row(suppmat, m), list));
od;
hull := convexhull(convert(pts, list));

for v from 1 to nops(hull) do
vertlist[v] := coordinates(hull[v]);
od;

eval(convert(vertlist, list));
end:
 
#input: list of vertices in counterclockwise order
#output: list of all lattice points in polygon
# "polygon fill" algorithm
# returns two lists: (points, interiorpoints);

getpoints := proc(hull)

local nvert, miny, minindex, maxy, maxindex, v, left,  leftstart,
rightstart, leftend, rightend, leftx, rightx, pts, intpts, numpts, numintpts, right, y, x, nextleft, nextright, leftslope, rightslope;


nvert := nops(hull);
miny := 32767;
maxy := -32767;
for v from 1 to nvert do
if ((hull[v])[2] < miny) then
   miny := (hull[v])[2];
   minindex := v;
fi;
if ((hull[v])[2] > maxy) then
   maxy := (hull[v])[2];
   maxindex := v;
fi;
od;


if  minindex = 1 then
  if (hull[nvert])[2] = miny then
      leftstart := nvert;
      rightstart := minindex;
  else
      if (hull[minindex+1])[2] = miny then
          leftstart := minindex;
          rightstart := minindex+1;
      else
      	leftstart := minindex;
      	rightstart := minindex;
      fi;
  fi;
else
  if (minindex = nvert) then
     leftstart := minindex;
     rightstart := minindex;
  else	
	if (hull[minindex+1])[2] = miny then
      		leftstart := minindex;	
      		rightstart := minindex+1;
	else
		leftstart := minindex;
      		rightstart := minindex;
	fi;
  fi;
fi;

v := rightstart;

while (hull[v][2] < maxy) do
  v := v + 1;
  if v > nvert then v := v - nvert; fi;
od;

rightend := v;

if (v < nvert) then
   if ((hull[v+1])[2] = maxy) then
      leftend := v+1;
   else
      leftend := v;
   fi;
else
   if ((hull[1])[2] = maxy) then
      leftend := 1;
   else
      leftend := v;
   fi;
fi;


v := leftstart;

left := leftstart;
nextleft := left - 1;
if (nextleft = 0) then nextleft := nvert; fi;
right := rightstart;
nextright := right + 1;
if (nextright > nvert) then nextright := 1; fi;
leftslope := (hull[nextleft][1] - hull[left][1])/(hull[nextleft][2] - hull[left][2]);
rightslope := (hull[nextright][1] - hull[right][1])/(hull[nextright][2] - hull[right][2]);

numpts := 0;
numintpts := 0;
leftx := hull[leftstart][1];
rightx := hull[rightstart][1];
for y from miny to maxy do
  for x from ceil(leftx) to floor(rightx) do
      numpts := numpts + 1;
      pts[numpts] := [x, y];
      
      if ((x > leftx) and (x < rightx) and (y > miny) and (y < maxy)) then
	  numintpts := numintpts + 1;
          intpts[numintpts] := [x, y];
      fi;
  od;
  if (y < maxy) then
  	if (y = hull[nextleft][2]) then
  	   left := nextleft;
     	   nextleft := left - 1;
           if  (nextleft = 0) then nextleft := nvert; fi;
     	   leftslope := (hull[nextleft][1] - hull[left][1])/(hull[nextleft][2] - hull[left][2]);
  	fi;
 	if (y = hull[nextright][2]) then
    	   right := nextright;
    	   nextright := right + 1;
           if  (nextright > nvert) then nextright := 1; fi;
           rightslope := (hull[nextright][1] - hull[right][1])/(hull[nextright][2] - hull[right][2]);
        fi;
        leftx := leftx + leftslope;
        rightx := rightx + rightslope;
  fi;            
od;
if (numintpts = 0) then
   intpts := []; fi;
return[convert(pts, list), convert(intpts, list)];

end:


### Given a list of vertices compute the normal fan and corresponding
### linear functional.

get_fan := proc(hull)

local nvert, fan, a, i, v, n, d;

nvert := nops(hull);
fan := array(1..nvert);
a := array(1..nvert);
for i from 1 to nvert-1 do
     v := (hull[i]-hull[i+1]);
     n := v[2];
     d := -v[1];
     n := n/gcd(n,d);
     d := d/gcd(v[2],d);
     fan[i] := [n,d];
     a[i] := -dotprod(vector(hull[i]), vector(fan[i]));
   od;
   v := (hull[nvert]-hull[1]);
   n := v[2]; 
   d := -v[1];
   n := n/gcd(n,d);
   d := d/gcd(v[2],d);
   fan[nvert] := [n,d];
   a[nvert] := -dotprod(vector(hull[nvert]), Vector(fan[nvert]));
   return(convert(eval(fan), list), convert(eval(a), list));
end:

get_homog := proc(pt, fan, a)

local va, v1;

v1 := map(v->dotprod(pt, vector(v)), fan);
v1 := v1 + a;

return(convert(v1, list));
end: 

add_supports := proc(supp1, supp2)

local size1, size2, i, j, S;

size1 := nops(supp1);
size2 := nops(supp2);
S := {};
for i from 1 to size1 do
    for j from 1 to size2 do
        S := S union {supp1[i] + supp2[j]};
    od;
od;
return convert(S, list);

end:

get_interior := proc(pts)

local numpts, size, interior, i, newpt, d;

interior := [];
numpts := nops(pts);
if (numpts = 0) then
  return interior
else
   size := nops(pts[1]);
   for i from 1 to size do
     d[i] := 1;
   od;
   d := convert(d, list);
   for i from 1 to numpts do
     if (verify(pts[i], d, 'list'('greater_equal'))) then
        interior := [op(interior), pts[i] - d];
     fi;
   od;
   return interior;
fi;
end:

Bezout_map := proc(R_1, R_2, R_3, alpha, supppts)

local Delta, e, w1, w, flag, i, v1, v, nmons, u1, u, beta, s, c;

Delta := 0; 
nmons := nops(supppts);
e := vector(nops(alpha), 1);
e := convert(e, list);
for w1 from 1 to nmons do
        w := supppts[w1];
        flag := false;
	for i in R_3 do
	if (i = 0) then
	   if ((w[1] + w[2]) >= (alpha[1] + alpha[2])) then
              flag := true;
              break;
           fi;
        else
           if (w[i] <= alpha[i]) then
               flag := true;
	       break;
	   fi;
        fi;
	od;
	    
	if (flag) then next; fi; 
	    for i in R_2 do
		if (w[i] <= alpha[i]) then
		    flag := true;
		    break;
		fi;
	    od;
	    if (not(flag)) then next; fi;
            
            for v1 from 1 to nmons do
 	    v := supppts[v1];
	    flag := false;
	    if (v1 = w1) then next; fi;
	    c := v + w;
	    for i in R_2 do
	        if (c[i] <= alpha[i]) then
		   flag := true;
		   break;
		fi;
	    od;
	    if (flag) then next; fi;
	    for i in R_1 do;
		if (c[i] <= alpha[i]) then 	
  	     	   flag := true;
		   break;
		fi;
	    od;
	    if (not(flag)) then next; fi;
           
	for u1 from 1 to nmons do	
            u := supppts[u1];
	    flag := false;
            if ((u1=v1) or (u1=w1)) then next; fi;
	    c := u + v + w;
	    for i in R_1 do;
		if (c[i] <= alpha[i]) then
		   flag := true;
		   break;
		fi;
	    od;
	    if (flag) then next; fi;
            
 
	    beta := c - alpha - e;
	    
	    
            #lprint(w1, v1, u1, alpha, beta);
            if (u1 <v1) then
	       if (v1 < w1) then
		  s := `[`||u1||`,`||v1||`,`||w1||`]`;
	       else
		  if (u1 < w1) then
	             s := -`[`||u1||`,`||w1||`,`||v1||`]`;
		  else
		     s := `[`||w1||`,`||u1||`,`||v1||`]`;
		  fi;
	       fi;
	     else
	       if (v1 > w1) then
	          s := -`[`||w1||`,`||v1||`,`||u1||`]`;
	       else
	         if (u1 < w1) then
	            s := -`[`||v1||`,`||u1||`,`||w1||`]`;
		 else
		    s := `[`||v1||`,`||w1||`,`||u1||`]`;
		 fi;
	       fi;
	     fi;      
	    
	    Delta := Delta + s*y^beta;
        od;
 od;
od;

return Delta;
end:

#wedgem := proc(Delta, qpts
bires := proc(supp, vars)

local suppmat, nmons, hull, fan, a, F, P, pts, intpts, qpts, intqpts, twoqpts, int2qpts, supppts ,inta, i, maxang, ang, eta_1, eta_2, M, R_1, R_2, R_3, c, B, e, m, j, n, alpha, w1, w, v1, v, u1, u, flag, beta, s;

supppts := [];

suppmat := supp_mat(supp, vars);
nmons := nops(supp);
for i from 1 to nmons do supppts := [op(supppts), row(suppmat, i)];
od;

hull := convhull(suppmat, nmons);
F := get_fan(hull);
fan := F[1];
a := F[2];
P := getpoints(hull);
pts := P[1];
intpts := P[2];

qpts := map(m->get_homog(m, fan, a), pts);
inta := map(x->x-1, a);
intqpts := map(m->get_homog(m, fan, inta), intpts);
twoqpts := add_supports(qpts, qpts);
int2qpts := get_interior(twoqpts);
supppts := map(m->get_homog(m, fan, a), supppts);

maxang := 0;

   #pick largest angle in fan
   for i from 1 to nops(fan)-1 do
	ang := angle(fan[i], fan[i+1]);
	if (evalf(ang) > evalf(maxang)) then
	   maxang := ang;
	   eta_1 := i;
	   eta_2 := i+1;
	fi;
    od;
    ang := angle(fan[nops(fan)], fan[1]);
    if (evalf(ang) > evalf(maxang)) then
       maxang := ang;
       eta_1 := nops(fan);
       eta_2 := 1;
    fi;
	
#partition fan

   M := inverse(matrix([fan[eta_1], fan[eta_2]]));
   R_1 := {};
   R_2 := {};
   R_3 := {};
  
   for i from 1 to nops(fan) do
       c := multiply(vector(fan[i]), M);
       if (c[1] >= 0 and c[2] <= 0) then
	  R_1 := R_1 union {i}
	else
	  if (c[2] >= 0) then
	     R_2 := R_2 union {i};
	  else
	     R_3 := R_3 union {i};
	  fi;
	fi;
   od;
  if (R_3 = {}) then R_3 := {0} fi;
    B := matrix(nops(int2qpts) + 3, nops(qpts) + 3*nops(intqpts), 0);
  e := vector(nops(fan), 1);
  e := convert(e, list);
  for i from 1 to nmons do
     m := supppts[i];
     member(m, qpts, 'cl');
     for j from 1 to 3 do	
        B[nops(int2qpts) + j, cl] := C[``||j||`,`||i];
        for n in intqpts do
            member(n+m, int2qpts, 'rw');
            member(n, intqpts, 'cl1');
            B[rw, nops(qpts) + (j-1)*nops(intqpts) + cl1] := C[``||j||`,`||i];
        od;
     od;
  od; 

  for alpha in qpts do
    member(alpha, qpts, 'cl');
    for w1 from 1 to nmons do
        w := supppts[w1];
        flag := false;
	for i in R_3 do
	if (i = 0) then
	   if ((w[1] + w[2]) >= (alpha[1] + alpha[2])) then
              flag := true;
              break;
           fi;
        else
           if (w[i] <= alpha[i]) then
               flag := true;
	       break;
	   fi;
        fi;
	od;
	    
	if (flag) then next; fi; 
	    for i in R_2 do
		if (w[i] <= alpha[i]) then
		    flag := true;
		    break;
		fi;
	    od;
	    if (not(flag)) then next; fi;
	
            for v1 from 1 to nmons do
 	    v := supppts[v1];
	    flag := false;
	    if (v1 = w1) then next; fi;
	    c := v + w;
	    for i in R_2 do
	        if (c[i] <= alpha[i]) then
		   flag := true;
		   break;
		fi;
	    od;
	    if (flag) then next; fi;
	    for i in R_1 do;
		if (c[i] <= alpha[i]) then 	
  	     	   flag := true;
		   break;
		fi;
	    od;
	    if (not(flag)) then next; fi;
	for u1 from 1 to nmons do	
            u := supppts[u1];
	    flag := false;
            if ((u1=v1) or (u1=w1)) then next; fi;
	    c := u + v + w;
	    for i in R_1 do;
		if (c[i] <= alpha[i]) then
		   flag := true;
		   break;
		fi;
	    od;
	    if (flag) then next; fi;
	    beta := c - alpha - e;
	    member(beta, int2qpts, 'rw');
	   
            #lprint(w1, v1, u1, alpha, beta);
            if (u1 <v1) then
	       if (v1 < w1) then
		  s := `[`||u1||`,`||v1||`,`||w1||`]`;
	       else
		  if (u1 < w1) then
	             s := -`[`||u1||`,`||w1||`,`||v1||`]`;
		  else
		     s := `[`||w1||`,`||u1||`,`||v1||`]`;
		  fi;
	       fi;
	     else
	       if (v1 > w1) then
	          s := -`[`||w1||`,`||v1||`,`||u1||`]`;
	       else
	         if (u1 < w1) then
	            s := -`[`||v1||`,`||u1||`,`||w1||`]`;
		 else
		    s := `[`||v1||`,`||w1||`,`||u1||`]`;
		 fi;
	       fi;
	     fi;      
	    
	    B[rw,cl] := B[rw,cl] + s;
        od;
 od;
od;
od;
eta_3 := eta_2 + 1;
  if (eta_3 > nops(fan)) then eta_3 := 1; fi;
 
 F_1 := {eta_1};
 F_2 := {eta_2};
 F_3 := {eta_3};
 if (eta_3 > eta_1) then 
 for i from (eta_3+1) to nops(fan) do 
     F_1 := F_1 union {i};
  od;
  for i from 1 to (eta_1 - 1) do
     F_1 := F_1 union {i};
  od;
 else
  for i from (eta_3+1) to (eta_1) do
    F_1 := F_1 union {i};
  od;
 fi;
  #lprint("R1 = ", R_1, "R2 = ", R_2, "R3 = ", R_3); 
#for alpha in twoqpts do
  #lprint("alpha");
  #Delta := Bezout_map(R_1, R_2, R_3, alpha, supppts);
  #lprint(Delta);
  #for beta in qpts do
     #member(beta, qpts, 'nm');
     #if member( alpha-beta, qpts) then
       #lprint(nm);
       #lprint(Bezout_map(R_1, R_2, R_3, alpha-beta, supppts));
    #fi;
 #od;
#od;
	
return evalm(B);
end:
