/* A cylinder.

  After loading the model, create the mesh with Modules->Mesh->3D, then in the options
  activate face coloring in Tools->Options->Mesh->Visibility->Surface faces
*/

Point(1) = {5,2,0,0.1};  // Center of the base of the cylinder

// 4 points to define the circle around Point 1
Point(2) = {5,2,-1,0.5};
Point(3) = {5,3,0,0.5};
Point(4) = {5,2,1,0.5};
Point(5) = {5,1,0,0.5};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
//Transfinite Surface{6} = {1,2,3,4};
Recombine Surface{6};
Physical Surface("bottom") = {6};

cyl[] = Extrude {-10,0,0} {
  Surface{6};
  Layers{100}; 
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

