/******************************************************************************
 A cylinder of radius 1 and length 2 along the y axis centered around {0,1,0}.
******************************************************************************/

ms = 0.4;                       // Mesh size

Point(1) = {0.5, 0, 0.5, 2*ms};  // Center of the base of the cylinder

// 4 points to define the circle around Point 1
Point(2) = { 1.5, 0, 0.5, ms};
Point(3) = { 0.5, 0, 1.5, ms};
Point(4) = {-0.5, 0, 0.5, ms};
Point(5) = { 0.5, 0, -0.5, ms};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
//Recombine Surface {6};

cyl[] = Extrude {0, 2, 0} {
  Surface{6};
  Layers { { 1,1,1,1 },
           { 1/20, 1/10, 1/5, 1}};
//  Recombine;  // Create hexahedral elements
};

//// Create groups for the definition of boundary conditions
// Remember that we have 4 sides (the Line Loop is made of 4 lines)

Physical Surface(1) = { 6 };        // base of the extruded entity
Physical Surface(2) = { cyl[0] };   // top of the extruded entity 
Physical Surface(3) = { cyl[2] };   // side 1
Physical Surface(4) = { cyl[3] };   // side 2
Physical Surface(5) = { cyl[4] };   // side 3
Physical Surface(6) = { cyl[5] };   // side 4
Physical Volume(1)  = { cyl[1] };   // extruded entity

