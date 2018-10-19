//Jamyd91's gmsh wrapper code version 0.2 (run the following in gmsh to generate a triangular mesh with topograpghy)
//2D mesh coordinates
cl=0.00;//define characteristic length
//Define surface points
Point(1) = {0.00,0.00,0.00,cl};//electrode location
//construct lines between each surface point
//add points below surface to make a fine mesh region
Point(2) = {0.00,0.00,0.00,cl*2.00};//base of smoothed mesh region
//make a polygon by defining lines between points just made.
//Add lines at the end points of each of the fine mesh region.
Line(1) = {1,2};
Line(2) = {1,2};
//compile lines into a line loop for a mesh surface/region.
Line Loop(1) = {1, -2};
Plane Surface(1) = {1};//Fine mesh region surface

//Add background region (Neumann boundary) points
cl2=0.00;//characteristic length for background region
Point(3) = {0.00,0.00,0.00,cl2};//far left upper point
Point(4) = {0.00,-0.00,0.00,cl2};//far left lower point
Point(5) = {0.00,0.00,0.00,cl2};//far right upper point
Point(6) = {0.00,-0.00,0.00,cl2};//far right lower point
//make lines encompassing all the background points - counter clock wise fashion
Line(3) = {1,3};
Line(4) = {3,4};
Line(5) = {4,6};
Line(6) = {6,5};
Line(7) = {5,1};
//Add line loops and plane surfaces to for nuemon region
Line Loop(2) = {3, 4, 5, 6, 7, 2, -1};
Plane Surface(2) = {2};//Fine mesh region surface

//Make a physical surface
Physical Surface(1) = {1, 2};

//Adding boreholes? 
//no more borehole strings to add.

//Adding buried electrodes? 
//no more buried electrode strings to add.

//Adding polygons?
//no more polygons to add.

//Adding boundaries?
//no more boundaries to add.

//j'ai fini!
