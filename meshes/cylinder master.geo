/* 
 A cylinder of radius 1 and length 10 along the z axis centered around {0,0,0}.
 The mesh is an extrusion with layers increasingly close as the middle of the
 cylinder is approached. The calculations for the Extrude command where done
 in a scheme session inside "README - cylinders.tm"
*/
Point(1) = {0,0,-5,0.5};

Point(2) = {1,0,-5,0.5};
Point(3) = {0,1,-5,0.5};
Point(4) = {-1,0,-5,0.5};
Point(5) = {0,-1,-5,0.5};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
Recombine Surface{6};
Physical Surface("bottom") = {6};

/*
  The Layers command was built using Scheme with a recursive definition of the
  sequence of points at which the layers where going to be. (see "README - cylinders.tm")
*/
cyl[] = Extrude {0,0,10} {
  Surface{6};
  Layers { { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 },
           {1/10, 7/38, 29/114, 605/1938, 116/323, 1922/4845, 413/969, 145/323, 301/646, 309/646, 1573/3230, 955/1938, 481/969, 161/323, 484/969, 2422/4845, 2423/4845, 485/969, 162/323, 488/969, 983/1938, 1657/3230, 337/646, 345/646, 178/323, 556/969, 2923/4845, 207/323, 1333/1938, 85/114, 31/38, 9/10, 1 }};
  Recombine;
};

Physical Surface("top") = {cyl[0]};
Physical Volume("cylinder") = {cyl[1]};
// Remeber that we have 4 sides (the Line Loop is made of 4 lines)
Physical Surface("side1") = {cyl[2]};
Physical Surface("side2") = {cyl[3]};
Physical Surface("side3") = {cyl[4]};
Physical Surface("side4") = {cyl[5]};

// Some presentation stuff
Color Red { Surface{6,cyl[0]}; }
Color Blue { Surface {cyl[2],cyl[3],cyl[4],cyl[5]};}
Geometry.Surfaces=1; // Won't work. Mabe because there's no 3D mesh yet?

// Store ids for use in files including this one (unused)
sidBottom = 6;
sidTop = cyl[0];
vidCylinder = cyl[1];
sidSide1 = cyl[2];
sidSide2 = cyl[3];
sidSide3 = cyl[4];
sidSide4 = cyl[5];

//Mesh.CharacteristicLengthFromCurvature = 1;

