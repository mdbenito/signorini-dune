/******************************************************************************
 A cylinder of radius 1 and length 10 along the x axis centered around {0,0,0}.
 The mesh is an extrusion with layers increasingly close as the middle of the
 cylinder is approached. The calculations for the Extrude command where done
 in a scheme session inside "README - cylinders.tm"
 UPD: I removed some of the middle layers of the extrusion. 
******************************************************************************/

ms = 0.4;                      // Mesh size

Point(1) = {5, 2, 0, 2*ms};  // Center of the base of the cylinder

// 4 points to define the circle around Point 1
Point(2) = {5, 2, -1, ms};
Point(3) = {5, 3,  0, ms};
Point(4) = {5, 2,  1, ms};
Point(5) = {5, 1,  0, ms};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
//Recombine Surface {6};

cylA[] = Extrude {-4, 0, 0} {
  Surface {6};
};

cylB[] = Extrude {-2, 0, 0} {
  Surface { cylA[0] };  // top of the first extrusion
  Layers { { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 },
           {1/10, 7/38, 29/114, 605/1938, 116/323, 1922/4845, 413/969, 145/323,
            301/646, 309/646, 1573/3230, 955/1938, 983/1938, 1657/3230,
            337/646, 345/646, 178/323, 556/969, 2923/4845, 207/323, 1333/1938,
            85/114, 31/38, 9/10, 1 }};
//  Recombine;  // Create hexahedral elements
};

cylC[] = Extrude { -4, 0, 0}{
  Surface { cylB[0] };   // top of the middle extrusion
};

//// Create groups for the definition of boundary conditions
// Remember that we have 4 sides (the Line Loop is made of 4 lines)

Physical Surface(1) = { 6 };         // base of the first extruded entity
Physical Surface(2) = { cylC[0] };   // top of the last extruded entity 
Physical Surface(3) = { cylA[4], cylA[5], cylC[4], cylC[5] }; // bottom side, no contact
Physical Surface(4) = { cylB[4], cylB[5] };                   // bottom side, contact
Physical Surface(5) = { cylA[2], cylA[3], cylB[2], cylB[3], cylC[2], cylC[3] };  // top side
Physical Volume(1)  = { cylA[1], cylB[1], cylC[1] };   // extruded entity

//// Some presentation stuff
// The colors seem to be ignored?

Color Red { Surface {6, cylC[0]}; }
Color Blue { Surface {cylB[2], cylB[3], cylB[2], cylB[3]};}

