//Gmsh wrapper code version 1.0 (run the following in gmsh to generate a triangular mesh with topograpghy)
//2D mesh coordinates
cl=0.11;//define characteristic length
//Define surface points
Point(1) = {-1.25,0.00,0.00,cl};//topography point
Point(2) = {0.00,0.00,0.00,cl};//electrode
Point(3) = {0.25,0.00,0.00,cl};//electrode
Point(4) = {0.50,0.00,0.00,cl};//electrode
Point(5) = {0.75,0.00,0.00,cl};//electrode
Point(6) = {1.00,0.00,0.00,cl};//electrode
Point(7) = {1.25,0.00,0.00,cl};//electrode
Point(8) = {1.50,0.00,0.00,cl};//electrode
Point(9) = {1.75,0.00,0.00,cl};//electrode
Point(10) = {2.00,0.00,0.00,cl};//electrode
Point(11) = {2.25,0.00,0.00,cl};//electrode
Point(12) = {2.50,0.00,0.00,cl};//electrode
Point(13) = {2.75,0.00,0.00,cl};//electrode
Point(14) = {3.00,0.00,0.00,cl};//electrode
Point(15) = {3.25,0.00,0.00,cl};//electrode
Point(16) = {3.50,0.00,0.00,cl};//electrode
Point(17) = {3.75,0.00,0.00,cl};//electrode
Point(18) = {4.00,0.00,0.00,cl};//electrode
Point(19) = {4.25,0.00,0.00,cl};//electrode
Point(20) = {4.50,0.00,0.00,cl};//electrode
Point(21) = {4.75,0.00,0.00,cl};//electrode
Point(22) = {5.00,0.00,0.00,cl};//electrode
Point(23) = {5.25,0.00,0.00,cl};//electrode
Point(24) = {5.50,0.00,0.00,cl};//electrode
Point(25) = {5.75,0.00,0.00,cl};//electrode
Point(26) = {7.00,0.00,0.00,cl};//topography point
//construct lines between each surface point
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13,14};
Line(14) = {14,15};
Line(15) = {15,16};
Line(16) = {16,17};
Line(17) = {17,18};
Line(18) = {18,19};
Line(19) = {19,20};
Line(20) = {20,21};
Line(21) = {21,22};
Line(22) = {22,23};
Line(23) = {23,24};
Line(24) = {24,25};
Line(25) = {25,26};
//add points below surface to make a fine mesh region
Point(27) = {-1.25,-2.50,0.00,cl*4.00};//base of smoothed mesh region
Point(28) = {0.25,-2.50,0.00,cl*4.00};//base of smoothed mesh region
Point(29) = {0.75,-2.50,0.00,cl*4.00};//base of smoothed mesh region
Point(30) = {1.25,-2.50,0.00,cl*4.00};//base of smoothed mesh region
Point(31) = {1.75,-2.50,0.00,cl*4.00};//base of smoothed mesh region
Point(32) = {2.25,-2.50,0.00,cl*4.00};//base of smoothed mesh region
Point(33) = {2.75,-2.50,0.00,cl*4.00};//base of smoothed mesh region
Point(34) = {3.25,-2.50,0.00,cl*4.00};//base of smoothed mesh region
Point(35) = {3.75,-2.50,0.00,cl*4.00};//base of smoothed mesh region
Point(36) = {4.25,-2.50,0.00,cl*4.00};//base of smoothed mesh region
Point(37) = {4.75,-2.50,0.00,cl*4.00};//base of smoothed mesh region
Point(38) = {5.25,-2.50,0.00,cl*4.00};//base of smoothed mesh region
Point(39) = {5.75,-2.50,0.00,cl*4.00};//base of smoothed mesh region
Point(40) = {7.00,-2.50,0.00,cl*4.00};//base of smoothed mesh region
//make a polygon by defining lines between points just made.
Line(26) = {27,28};
Line(27) = {28,29};
Line(28) = {29,30};
Line(29) = {30,31};
Line(30) = {31,32};
Line(31) = {32,33};
Line(32) = {33,34};
Line(33) = {34,35};
Line(34) = {35,36};
Line(35) = {36,37};
Line(36) = {37,38};
Line(37) = {38,39};
Line(38) = {39,40};
//Add lines at the end points of each of the fine mesh region.
Line(39) = {1,27};
Line(40) = {26,40};
//compile lines into a line loop for a mesh surface/region.
Line Loop(1) = {39, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, -40, -25, -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1};
Plane Surface(1) = {1};//Fine mesh region surface

//Add background region (Neumann boundary) points
cl2=5.47;//characteristic length for background region
Point(41) = {-26.25,0.00,0.00,cl2};//far left upper point
Point(42) = {-26.25,-27.50,0.00,cl2};//far left lower point
Point(43) = {32.00,0.00,0.00,cl2};//far right upper point
Point(44) = {32.00,-27.50,0.00,cl2};//far right lower point
//make lines encompassing all the background points - counter clock wise fashion
Line(41) = {1,41};
Line(42) = {41,42};
Line(43) = {42,44};
Line(44) = {44,43};
Line(45) = {43,26};
//Add line loops and plane surfaces to for nuemon region
Line Loop(2) = {41, 42, 43, 44, 45, 40, -38, -37, -36, -35, -34, -33, -32, -31, -30, -29, -28, -27, -26, -39};
Plane Surface(2) = {2};//Fine mesh region surface

//Make a physical surface
Physical Surface(1) = {1, 2};

//Adding polygons?
//end of polygons.

//Adding boundaries?
//end of boundaries.

//j'ai fini!
