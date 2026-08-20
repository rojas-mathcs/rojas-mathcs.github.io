# Copyright 2011, Jeff Sommars.
# This code represents part of a final project for J. Maurice Rojas' 
# summer 2011 REU class. Mounir Nisse, J. Maurice Rojas, and Korben Rusek 
# co-advised Sommars in summer 2011. Many thanks to Daniel Smith. 

def vectorcmp(X,Y):
#Used to sort by slope for sortedmatrix. Necessary for removing vectors with duplicate slopes.
	x = slope(X)
	y = slope(Y)
	if x < y:
		return int(-1)
	if x > y:
		return int(1)
	return int(0)


def upsortcmp(X,Y):
#Used to sort by trianglepointfinder when making the adjacency matrix. Sorts smallest y value to largest y value.
	if X[1] < Y[1]:
		return int(-1)
	elif X[1] > Y[1]:
		return int(1)
	return int(0)


def rightsortcmp(X,Y):
#Used to sort by trianglepointfinder when making the adjacency matrix. Sorts smallest x value to largest x value.
	if X[0] < Y[0]:
		return int(1)
	if X[0] > Y[0]:
		return int(-1)
	return int(0)


def slope(X):
#Finds the slope of a vector, called by vectorcmp.
	x = X[0]
	y = X[1]
	if y == 0:
		return -oo
	slope = -x/y
	return slope


def sortedmatrix(A):
#First, it finds the bmatrix and sorts it as described in the algorithm. It then combines vectors in the bmatrix with the same slope. These vectors could cause trouble later on in the program, so it's prudent to remove them early.
	a = len(A[0])
	B = A.right_kernel().echelonized_basis_matrix().transpose()
	C = list(B)
	C.sort(vectorcmp, reverse=true)
	for k in range(0,a-1):
		for j in range(1,a):
			if (k<j):
				if (C[k][1] == 0 and C[j][1] == 0 and C[k][0] != 0 and C[j][0] != 0):
					C[k] = C[k] + C[j]
					C[j] = (0,0)
				if (C[k][1] != 0 and C[j][1] != 0):
					if (C[k][0]/C[k][1] == C[j][0]/C[j][1]):
						C[k] = C[k] + C[j]
						C[j] = (0,0)
	goodbmatrix = []
	bmatrixlength = 0
	for k in range(0,a):
		if(C[k][0] != 0 or C[k][1] != 0):
			goodbmatrix += [C[k]]
			bmatrixlength += 1
	return goodbmatrix, bmatrixlength


def kapranov(A):
#Uses the Horn-Kapranov parametrization to find the starting point of the coamoeba.
	B = sortedmatrix(A)[0]
	N = sortedmatrix(A)[1]
	var('t')
	first = 1
	second = 1
	for k in range(0, N):
		first = first*((B[k][0]+t*B[k][1]))^(B[k][0])
		second = second*((B[k][0]+t*B[k][1]))^(B[k][1])
	x = -1
	y = -1
	first = first > 0
	second = second > 0
	j = sortedmatrix(A)[1]
	if (j == 0):
		print "Sorry, your bmatrix is empty, so there is no coamoeba to draw."
		return
	k = sortedmatrix(A)[0][0][0]
	m = sortedmatrix(A)[0][0][1]
	big = k/m + 1
	if bool(first.subs(t == big)):
		x = 1
	if bool(second.subs(t == big)):
		y = 1

	#converts the x and y values to the Po value where the coamoeba starts	
	M = MatrixSpace(RR,1,2)
	if x == 1 and y == 1:
		M = [[0,0]]
	elif x == 1 and y == -1:
		M = [[0,pi]]
	elif x == -1 and y == 1:
		M = [[pi,0]]
	elif x == -1 and y == -1:
		M = [[pi,pi]]
	return M


def maxdimensions(A):
#Finds the greatest magnitude of x and y that the principle coamoeba reaches.
	bigx = 0
	bigy = 0
	l = sortedmatrix(A)[1]
	if (l == 0):
		print "Sorry, your bmatrix is empty, so there is no coamoeba to draw."
		return
	points = MatrixSpace(RR,l+1,2)(range((l+1)*2))
	points[0] = kapranov(A)[0]
	bmatrix = sortedmatrix(A)[0]
	for k in range(1,l+1):
		points[k] = points[k-1] + pi*bmatrix[k-1]
	for k in range(0,l+1):
		i = abs(points[k][0])
		j = abs(points[k][1])
		if(bigx < i):
			bigx = i
		if(bigy < j):
			bigy = j
	bigmatrix = MatrixSpace(RR,1,2)
	bigmatrix = [[bigx,bigy]]
	return bigmatrix


def xrange(A):
#Gives all the x coordinates where the prinicple coamoeba should be drawn in order to get the fundamental domain correct.
	xstart = -maxdimensions(A)[0][0]
	xstart = xstart - (float(xstart)%float(2*pi))
	xinc = float(abs(xstart)*2)/float(2*pi)+1
	xinc = int(xinc)
	xvalues = MatrixSpace(RR,xinc+2,1)(range(xinc+2))
	xvalues[0] = [xstart - 2*pi]
	incmatrix = MatrixSpace(RR,1,1)([[2*pi]])
	for k in range(1,xinc+2):
		xvalues[k] = xvalues[k-1] + incmatrix[0]
	return xvalues


def yrange(A):
#Gives all the y coordinates where the prinicple coamoeba should be drawn in order to get the fundamental domain correct.
	ystart = -maxdimensions(A)[0][1]
	ystart = ystart - (float(ystart)%float(2*pi))
	yinc = float(abs(ystart)*2)/float(2*pi)+1
	yinc = int(yinc)
	yvalues = MatrixSpace(RR,yinc+2,1)(range(yinc+2))
	yvalues[0] = [ystart - 2*pi]
	incmatrix = MatrixSpace(RR,1,1)([[2*pi]])
	for k in range(1,yinc+2):
		yvalues[k] = yvalues[k-1] + incmatrix[0]
	return yvalues



def coamoeba(A, centered=true, principle=false, minx=0, maxx=0, miny=0, maxy=0, mult=false, x=0, y=0):
	"""
    	   This is the function that should be called for drawing two dimensional coamoebas.
	   It takes an A matrix as an input, which should be given to Sage as a MatrixSpace with values in Z.

	   By default, this function will return a graph of the fundamental domain -pi<=x,y<=pi.
	   If centered=false, the graph will be shifted so that its center is at (pi,pi).
	   If principle=true, a single coamoeba will be drawn that will not stay in the fundamental domain. This is particularly good for seeing the distinct chambers of the coamoeba.
	   If centered=true, setting minx, maxx, miny and maxy such that -pi<=minx<maxx<=pi and -pi<=miny<maxy<=pi allows you to view a smaller area of the fundamental domain -pi<=x,y<=pi in higher detail.
	   Finally, you can set mult=true and pass in x and y values to receive the multiplicity at that point, as long as it is within -pi<x<pi and -pi<y<pi. Do not test the multiplicity on any boundary of the coamoeba, interior or exterior.
	   
	   For testing purposes, a few matrices have been included. They are:

A = MatrixSpace(ZZ,1,3)([[1,1,1]])

B = MatrixSpace(ZZ,2,4)([[1,1,1,1],[0,1,2,3]])

C = MatrixSpace(ZZ,3,5)([[1,1,1,1,1],[2,1,0,2,3],[3,2,0,1,2]])

D = MatrixSpace(ZZ,4,6)([[1,1,1,1,1,1],[0,1,0,1,0,1],[0,0,1,1,0,0],[0,0,0,0,1,1]])

E = MatrixSpace(ZZ,6,8)([[1,1,1,1,1,1,1,1],[-2,-3,0,1,0,0,0,0],[3,4,0,0,1,0,0,0],[-4,-5,0,0,0,1,0,0],[14,16,0,0,0,0,1,0],[-11,-13,0,0,0,0,0,1]])

	   A, B, and C are the examples in the paper by Nilsson and Passare, D is a degenerate case with repeated roots, and E is an interesting extremal case. Allow E some time to run.

	   The primary algorithm implemented can be found in the paper "Discriminanat coamoebas in dimension two" by Lisa Nilsson and Mikael Passare.
	   Much thanks to Dr. Rojas, Dr. Nisse, Korben Rusek, and Daniel Smith; each of them helped me greatly.

	   This program was written by Jeff Sommars (Wheaton College) as part of the 2011 Math REU at Texas A&M
	"""


	l = sortedmatrix(A)[1]
	if (l == 0):
		print "Sorry, your bmatrix is empty, so there is no coamoeba to draw."
		return
	X = xrange(A)
	Y = yrange(A)
	i = len(X.transpose()[0])
	j = len(Y.transpose()[0])

	points = MatrixSpace(RR,2*l+1,2)(range((2*l+1)*2))
	points[0] = kapranov(A)[0]
	bmatrix = sortedmatrix(A)[0]
	for k in range(1,l+1):
		points[k] = points[k-1] + pi*bmatrix[k-1]
	for k in range(l+1,2*l+1):
		points[k] = points[k-1] - pi*bmatrix[k-1-l]


	#Gives the multiplicity of a point in the fundamental domain from 0 to 2pi or equivalently from -pi to pi.
	if(mult==true):
		pointstester = MatrixSpace(RR,2*l+1,2)(range(4*l+2))
		multipl = 0
		for k in range(0,i):
			for z in range(0,j):
				for q in range(0,2*l+1):
					pointstester[q] = (points[q][0]+X[k][0],points[q][1]+Y[z][0])
				multipl += multiplicity(pointstester,x,y)
		return multipl



	testpoints = trianglepointfinder(points,l)

	#This draws the principle coamoeba alone, if that is what is desired.
	if(principle==true):
		primarygraph = polygon2d([0,0])
		for k in range(0,len(testpoints)):
			pointstester = copy(points)
			multipl = multiplicity(pointstester,testpoints[k][6],testpoints[k][7])
			primarygraph+=polygon2d(([testpoints[k][0],testpoints[k][1]],[testpoints[k][2],testpoints[k][3]],[testpoints[k][4],testpoints[k][5]]),alpha=min(.15*multipl,1))
		show(primarygraph,aspect_ratio=1)
		return



	#This finds all the coamoebas that intersect the fundamental domain centered at (0,0) and removes all the coamoebas that do not intersect the fundamental domain.
	if(centered==true):
		matrixlist = []
		incmat = MatrixSpace(RR,2,1)(range(2))
		points2 = MatrixSpace(RR,2*l+1,2)(range((2*l+1)*2))
		for k in range(0,i):
			for p in range(0,j):
				incmat[0] = X[k]
				incmat[1] = Y[p]
				incmat = incmat.transpose()
				for m in range (0,2*l+1):
					points2[m] = points[m] + incmat[0]
				for c in range(0,2*l):
					x1 = points2[c][0]
					y1 = points2[c][1]
					x2 = points2[c+1][0]
					y2 = points2[c+1][1]
					slope = (y2-y1)/(x2-x1)
					if(-pi<=x1<=pi and -pi<=y1<=pi):
						incmat = incmat.transpose()
						matrixlist.append((incmat[0][0],incmat[1][0]))
						break
					elif(slope*(-3.1415926535979323-x1)+y1>-pi and slope*(-3.1415926535979323-x1)+y1<pi and ((x1<-pi and x2>-pi) or (x1>-pi and x2<-pi))):
						incmat = incmat.transpose()
						matrixlist.append((incmat[0][0],incmat[1][0]))
						break
					elif(slope*(3.1415926535979323-x1)+y1>-pi and slope*(3.1415926535979323-x1)+y1<pi  and ((x1<pi and x2>pi) or (x1>pi and x2<pi))):
						incmat = incmat.transpose()
						matrixlist.append((incmat[0][0],incmat[1][0]))
						break
					elif((3.1415926535979323-y1)/slope+x1>-pi and (3.1415926535979323-y1)/slope+x1<pi and ((y1<-pi and y2>-pi) or (y1>-pi and y2<-pi))):
						incmat = incmat.transpose()
						matrixlist.append((incmat[0][0],incmat[1][0]))
						break
					elif((-3.1415926535979323-y1)/slope+x1>-pi and (-3.1415926535979323-y1)/slope+x1<pi)  and ((y1<pi and y2>pi) or (y1>pi and y2<pi)):
						incmat = incmat.transpose()
						matrixlist.append((incmat[0][0],incmat[1][0]))
						break
					elif(c ==(2*l-1)):
						incmat = incmat.transpose()



	#This finds all the coamoebas that intersect the fundamental domain centered at (pi,pi) and removes all the coamoebas that do not intersect the fundamental domain.
	if(centered==false):
		matrixlist = []
		incmat = MatrixSpace(RR,2,1)(range(2))
		points2 = MatrixSpace(RR,2*l+1,2)(range((2*l+1)*2))
		for k in range(0,i):
			for p in range(0,j):
				incmat[0] = X[k]
				incmat[1] = Y[p]
				incmat = incmat.transpose()
				for m in range (0,2*l+1):
					points2[m] = points[m] + incmat[0]
				for c in range(0,2*l):
					x1 = points2[c][0]
					y1 = points2[c][1]
					x2 = points2[c+1][0]
					y2 = points2[c+1][1]
					slope = (y2-y1)/(x2-x1)
					if(0<=x1<=2*pi and 0<=y1<=2*pi):
						incmat = incmat.transpose()
						matrixlist.append((incmat[0][0],incmat[1][0]))
						break
					elif(slope*(-x1)+y1>0 and slope*(-x1)+y1<2*pi and ((x1<-pi and x2>-pi) or (x1>-pi and x2<-pi))):
						incmat = incmat.transpose()
						matrixlist.append((incmat[0][0],incmat[1][0]))
						break
					elif(slope*(2*3.1415926535979323-x1)+y1>0 and slope*(2*3.1415926535979323-x1)+y1<2*pi  and ((x1<pi and x2>pi) or (x1>pi and x2<pi))):
						incmat = incmat.transpose()
						matrixlist.append((incmat[0][0],incmat[1][0]))
						break
					elif((2*3.1415926535979323-y1)/slope+x1>0 and (2*3.1415926535979323-y1)/slope+x1<2*pi and ((y1<-pi and y2>-pi) or (y1>-pi and y2<-pi))):
						incmat = incmat.transpose()
						matrixlist.append((incmat[0][0],incmat[1][0]))
						break
					elif((-y1)/slope+x1>0 and (-y1)/slope+x1<2*pi)  and ((y1<pi and y2>pi) or (y1>pi and y2<pi)):
						incmat = incmat.transpose()
						matrixlist.append((incmat[0][0],incmat[1][0]))
						break
					elif(c ==(2*l-1)):
						incmat = incmat.transpose()



	#This tests the points in the triangles and adds the graphs of each individual triangle (colored appropriately) to the final graph.
	finalgraph = polygon2d([0,0])
	initialpointslist=list(points)
	sendmatrix = MatrixSpace(RR,2*l+1,2)(range(4*l+2))
	for k in range(0,len(matrixlist)):
		for p in range(0,2*l+1):
			sendmatrix[p]=(initialpointslist[p][0]+matrixlist[k][0],initialpointslist[p][1]+matrixlist[k][1])
		for z in range(0,len(testpoints)):
			x1=testpoints[z][0]+matrixlist[k][0]
			y1=testpoints[z][1]+matrixlist[k][1]
			x2=testpoints[z][2]+matrixlist[k][0]
			y2=testpoints[z][3]+matrixlist[k][1]
			x3=testpoints[z][4]+matrixlist[k][0]
			y3=testpoints[z][5]+matrixlist[k][1]
			x4=testpoints[z][6]+matrixlist[k][0]
			y4=testpoints[z][7]+matrixlist[k][1]
			multipl = multiplicity(copy(sendmatrix),x4,y4)
			if (multipl!=0):
				finalgraph += (polygon2d(([x1,y1],[x2,y2],[x3,y3]),alpha=min(multipl*.15,1)))


	if(centered==true):
		if(minx!=0 or maxx!=0 or miny!=0 or maxy!=0):
			show(finalgraph, xmin=minx, xmax=maxx, ymin=miny, ymax=maxy, aspect_ratio=1)
			return
		else:
			show(finalgraph, xmin = -3.021, xmax = 3.021, ymin = -3.021, ymax = 3.021, aspect_ratio=1)
			return
	if(centered==false):
		show(finalgraph, xmin = 0.12, xmax = 2*pi-.12, ymin = 0.12, ymax = 2*pi-.12, aspect_ratio=1)
		return




def trianglepointfinder(points,l):
#This algorithm takes the set of points that are obtained from the algorithm by Passare and Nilsson. This set of points defines a complex polygon. trianglepointfinder takes the convex hull of the complex polygon, and subdivides it into triangles such that no side of a triangle intersects one of the complex polygons sides or the side of another triangle. It returns the triangles and points within each triangle. This is done so that the point in each triangle can be tested to find the multiplicity of the triangle, and then the triangle can be drawn with the correct multiplicity. Note that there is some small possibility of rounding error wherever the round function is used.


	#It starts by finding all the points where any line segments intersect each other
	intersectionlist = []
	for k in range(0,2*l):
		for j in range(1,2*l):
			if(k<(j) and j!=k+1):
				x1 = points[k][0]
				y1 = points[k][1]
				x2 = points[k+1][0]
				y2 = points[k+1][1]
				x3 = points[j][0]
				y3 = points[j][1]
				x4 = points[j+1][0]
				y4 = points[j+1][1]
				a = ((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
				b = ((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
				if(min(x1,x2)<=a<=max(x1,x2) and min(x3,x4)<=a<=max(x3,x4) and min(y1,y2)<=b<=max(y1,y2) and min(y3,y4)<=b<=max(y3,y4) and not (math.isnan(a) or math.isnan(b))):
					intersectionlist += [(a,b)]


	#It removes all the points of intersection
	intersectionlistwithoutduplicates = list(set(copy(intersectionlist)))
	intersectionlistwithoutduplicates.sort(upsortcmp)
	intersectionlistsize = len(intersectionlistwithoutduplicates)


	#This makes fullpointlist, a list that contains all the points of intersection and all of points, without any duplicates.
	fullpointlistwithduplicates = list(points)
	fullpointlistwithduplicates += intersectionlistwithoutduplicates
	fullpointlistwithduplicates.sort(upsortcmp)
	for k in range(0,2*l+intersectionlistsize):
		for j in range(0,2*l+intersectionlistsize+1):
			if(k<j):
				if(fullpointlistwithduplicates[k][0]==fullpointlistwithduplicates[j][0] and fullpointlistwithduplicates[k][1]==fullpointlistwithduplicates[j][1]):
					fullpointlistwithduplicates[j] = (infinity,0)
	fullpointlist = []
	for k in range(0,2*l+intersectionlistsize+1):
		if(fullpointlistwithduplicates[k][0]!=infinity and fullpointlistwithduplicates[k][0]!=2.0/0.0):
			fullpointlist += [fullpointlistwithduplicates[k]]


	#These are used when creating the adjacency matrix.
	up = copy(intersectionlistwithoutduplicates)
	intersectionlistwithoutduplicates.sort(upsortcmp,reverse=true)
	down = copy(intersectionlistwithoutduplicates)
	intersectionlistwithoutduplicates.sort(rightsortcmp)
	left = copy(intersectionlistwithoutduplicates)
	intersectionlistwithoutduplicates.sort(rightsortcmp,reverse=true)
	right = copy(intersectionlistwithoutduplicates)


	#From here, I make an adjacency matrix. I start by making the matrix and then I fill in "1's" where appropriate.
	adjacencymatrix = []
	adjacrow = []
	for k in range(0,len(fullpointlist)):
		adjacrow += [0]
	for k in range(0,len(fullpointlist)):
		adjacencymatrix += [copy(adjacrow)]
	
	tester = []
	for i in range(len(fullpointlist)):
		tester += [list(fullpointlist[i])]
	for k in range(0,2*l):
		x1 = points[k][0]
		y1 = points[k][1]
		x2 = points[k+1][0]
		y2 = points[k+1][1]
		slope=(y2-y1)/(x2-x1)
		vert=x2-x1
		if(y1<y2):
			for j in range(0,intersectionlistsize+1):
				if(j==intersectionlistsize):
					xcoord = tester.index([x1,y1])
					ycoord = tester.index([x2,y2])
					adjacencymatrix[xcoord][ycoord] = 1
					break
				if((up[j][1]>y1 and up[j][1]<y2 and round(up[j][1]-y2,13)==round(slope*(up[j][0]-x2),13)) or (vert==0 and up[j][1]>y1 and up[j][1]<y2)):
					xcoord = tester.index([up[j][0],up[j][1]])
					ycoord = tester.index([x1,y1])
					x1=up[j][0]
					y1=up[j][1]
					adjacencymatrix[xcoord][ycoord] = 1
				elif(up[j][1]>=y2):
					xcoord = tester.index([x1,y1])
					ycoord = tester.index([x2,y2])
					adjacencymatrix[xcoord][ycoord] = 1
					break
		elif(y1>y2):
			for j in range(0,intersectionlistsize+1):
				if(j==intersectionlistsize):
					xcoord = tester.index([x1,y1])
					ycoord = tester.index([x2,y2])
					adjacencymatrix[xcoord][ycoord] = 1
					break
				if((down[j][1]<y1 and down[j][1]>y2 and round(down[j][1]-y2,13)==round(slope*(down[j][0]-x2),13)) or (vert==0 and down[j][1]<y1 and down[j][1]>y2)):
					xcoord = tester.index([down[j][0],down[j][1]])
					ycoord = tester.index([x1,y1])
					x1=down[j][0]
					y1=down[j][1]
					adjacencymatrix[xcoord][ycoord] = 1
				elif(down[j][1]<=y2):
					xcoord = tester.index([x1,y1])
					ycoord = tester.index([x2,y2])
					adjacencymatrix[xcoord][ycoord] = 1
					break
		elif(x1<x2):
			for j in range(0,intersectionlistsize+1):
				if(j==intersectionlistsize):
					xcoord = tester.index([x1,y1])
					ycoord = tester.index([x2,y2])
					adjacencymatrix[xcoord][ycoord] = 1
					break
				if(right[j][0]>x1 and right[j][0]<x2 and round(right[j][1]-y2,13)==round(slope*(right[j][0]-x2),13)):
					xcoord = tester.index([right[j][0],right[j][1]])
					ycoord = tester.index([x1,y1])
					x1=right[j][0]
					y1=right[j][1]
					adjacencymatrix[xcoord][ycoord] = 1
				elif(right[j][0]>x2):
					xcoord = tester.index([x1,y1])
					ycoord = tester.index([x2,y2])
					adjacencymatrix[xcoord][ycoord] = 1
					break
		elif(x1>x2):
			for j in range(0,intersectionlistsize+1):
				if(j==intersectionlistsize):
					xcoord = tester.index([x1,y1])
					ycoord = tester.index([x2,y2])
					adjacencymatrix[xcoord][ycoord] = 1
					break
				if(left[j][0]<x1 and left[j][0]>x2 and round(left[j][1]-y2,13)==round(slope*(left[j][0]-x2),13)):
					xcoord = tester.index([left[j][0],left[j][1]])
					ycoord = tester.index([x1,y1])
					x1=left[j][0]
					y1=left[j][1]
					adjacencymatrix[xcoord][ycoord] = 1
				elif(left[j][0]<x2):
					xcoord = tester.index([x1,y1])
					ycoord = tester.index([x2,y2])
					adjacencymatrix[xcoord][ycoord] = 1
					break


	#This makes the adjacency matrix uppertriangular.
	for k in range(0,len(fullpointlist)):
		for j in range(0,len(fullpointlist)):
			if(adjacencymatrix[k][j]!=0 and k>j):
				adjacencymatrix[k][j]=0
				adjacencymatrix[j][k]=1


	#This determines where all the lines actually exist. The indices correspond to fullpointlist.
	hasline = []
	for k in range(0,len(fullpointlist)-1):
		for j in range(0,len(fullpointlist)):
			if(j>k and adjacencymatrix[k][j]!=0):
				hasline += [[k,j]]


	#This triangulates the graph of the coamoeba and adds to the adjacency matrix 1s for all the lines that are created.
	for k in range(0,len(fullpointlist)):
		for j in range(0,len(fullpointlist)):
			if(adjacencymatrix[k][j]==0 and j>k):
				for m in range(0,len(hasline)):
					x1 = fullpointlist[k][0]
					y1 = fullpointlist[k][1]
					x2 = fullpointlist[j][0]
					y2 = fullpointlist[j][1]
					x3 = fullpointlist[hasline[m][0]][0]
					y3 = fullpointlist[hasline[m][0]][1]
					x4 = fullpointlist[hasline[m][1]][0]
					y4 = fullpointlist[hasline[m][1]][1]
					slope1 = (y2-y1)/(x2-x1)
					slope3 = (y4-y3)/(x4-x3)
					s = ((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))/((x4-x3)*(y2-y1)-(y4-y3)*(x2-x1))
					t = (x3-x1+(x4-x3)*s)/(x2-x1)
					if (0<round(s,13)<1 and 0<round(t,13)<1):
						break
					if(round((y1-slope1*(x1)),13)==round((y3-slope3*(x3)),13) and round(slope1,13)==round(slope3,13) and abs(slope1)!=oo):
						if(x1<=x3<=x2 or x1<=x4<=x2 or x2<=x3<=x1 or x2<=x4<=x1):
							break
					if(abs(slope1)==abs(slope3) and abs(slope1)==oo):
						if(round(x1,13)==round(x3,13)):
							if(y1<=y3<=y2 or y1<=y4<=y2 or y2<=y3<=y1 or y2<=y4<=y1):
								break
					if(m == len(hasline)-1 and k!=j):
						hasline += [[k,j]]
						adjacencymatrix[k][j] += 1
						break


	#This simple algorithm searches through the adjacency matrix to figure out where there are triangles formed. For a triangle to exist, it is necessary for (i,j) or (j,i) to equal 1, (i,k) or (k,i) to equal 1, and (k,j) or (j,k) to equal 1. This means that there would be lines connecting each of points i j and k to each other.
	trianglelist = []
	for k in range(0,len(fullpointlist)):
		for j in range(0,len(fullpointlist)):
			for i in range(0,len(fullpointlist)):
				if(j>k and adjacencymatrix[k][j]==1 and i>j and adjacencymatrix[k][i]==1):
					if(adjacencymatrix[j][i]==1 or adjacencymatrix[i][j]==1):
						trianglelist += [[i,j,k]]

	


	#Once the triangles have been found, a point inside each triangle is found so that it can be sent to the multiplicity checker. Once a point is found (x4,y4) that is contained within the triangle, each of the vertices along with the interior point are made into a small list and added onto sendlist, the return of the function.
	sendlist = []
	for k in range(0,len(trianglelist)):
		x1=fullpointlist[trianglelist[k][0]][0]
		y1=fullpointlist[trianglelist[k][0]][1]
		x2=fullpointlist[trianglelist[k][1]][0]
		y2=fullpointlist[trianglelist[k][1]][1]
		x3=fullpointlist[trianglelist[k][2]][0]
		y3=fullpointlist[trianglelist[k][2]][1]
		x4=(x3+(x1+x2)/2)/2
		y4=(y3+(y1+y2)/2)/2
		sendlist += [[x1,y1,x2,y2,x3,y3,x4,y4]]

#	goodtester(adjacencymatrix,fullpointlist,len(fullpointlist))

	return sendlist

# def goodtester(A,goodp,jj):
#    grapher = line(([0,0],[0,0]))
#    for k in range(0,jj):
#        for j in range(0,jj):
#            if A[k][j]==1:
#                grapher+=line((goodp[k],goodp[j]))
#    show(grapher)
#    return

def multiplicity(A,x,y):
#This algorithm for finding the multiplicity of a point was implemented from the paper "A Winding Number and Point-in-Polygon Algorithm" by Rick Miranda of Colorado State University.
#As a warning, it may not return the correct multiplicity if the point chosen is on an edge, interior or exterior. It may return an answer one less than the true muliplicity for every edge on which the point lies.
	points = A
	l = len(points.transpose()[0])
	shift = MatrixSpace(RR,1,2)([[x,y]])
	for k in range(0,l):
		points[k] = points[k]-shift[0]
	mult = 0
	for k in range(0,l-1):
		if(points[k][1]*points[k+1][1] < 0):
			r = points[k][0] + points[k][1]*(points[k+1][0]-points[k][0])/(points[k][1]-points[k+1][1])
			if(r>0):
				if(points[k][1]<0):
					mult = mult + 1
				else:
					mult = mult - 1
		elif(points[k][1] == 0 and points[k][0] >= 0):
			if(points[k+1][1]>0):
				mult = mult+1/2
			else:
				mult = mult-1/2
		elif(points[k+1][1] == 0 and points[k+1][0] >= 0):
			if(points[k][1] < 0):
				mult = mult + 1/2
			else:
				mult = mult - 1/2
	if(mult < 0):
		mult = 0
	return mult
