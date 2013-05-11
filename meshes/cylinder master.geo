/******************************************************************************
 A cylinder of radius 1 and length 10 along the z axis centered around {0,0,0}.
 The mesh is an extrusion with layers increasingly close as the middle of the
 cylinder is approached. The calculations for the Extrude command where done
 in a scheme session inside "README - cylinders.tm"
******************************************************************************/

ms = 0.8;                       // Mesh size

Point(1) = {0, 0, -5, 2*ms};  // Center of the base of the cylinder

// 4 points to define the circle around Point 1
Point(2) = { 1,  0, -5, ms};
Point(3) = { 0,  1, -5, ms};
Point(4) = {-1,  0, -5, ms};
Point(5) = { 0, -1, -5, ms};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
//Recombine Surface {6};

cyl[] = Extrude {0, 0, 10} {
  Surface{6};
  Layers { { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 },
           {1/10, 7/38, 29/114, 605/1938, 116/323, 1922/4845, 413/969, 145/323,
            301/646, 309/646, 1573/3230, 955/1938, 481/969, 
            488/969, 983/1938, 1657/3230,
            337/646, 345/646, 178/323, 556/969, 2923/4845, 207/323, 1333/1938,
            85/114, 31/38, 9/10, 1 }};
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

//// Some presentation stuff
// The colors seem to be ignored?

Color Red { Surface {6, cyl[0]}; }
Color Blue { Surface {cyl[2], cyl[3], cyl[4], cyl[5]};}

//// Store ids for use in files including this one (unused)

sidBottom = 6;
sidTop = cyl[0];
vidCylinder = cyl[1];
sidSide1 = cyl[2];
sidSide2 = cyl[3];
sidSide3 = cyl[4];
sidSide4 = cyl[5];

