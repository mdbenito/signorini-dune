/* A prism extruded from a square */

ms = 0.4;     // Mesh size
height = 1.0; // Prism height

Point(1) = {0, -0.01, 0, ms};
Point(2) = {1, -0.01, 0, ms};
Point(3) = {1, -0.01, 1, ms};
Point(4) = {0, -0.01, 1, ms};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
//Transfinite Surface {6};
//Recombine Surface {6};

cube[] = Extrude {0, -height, 0} { 
  Surface {6};
  Layers {10*(1-ms)};
  //Recombine;
};

Physical Volume(1)  = {cube[1]}; // cube
Physical Surface(1) = {6};       // base
Physical Surface(2) = {cube[2]}; // side
Physical Surface(3) = {cube[3]}; // side
Physical Surface(4) = {cube[4]}; // side
Physical Surface(5) = {cube[5]}; // side
Physical Surface(6) = {cube[0]}; // bottom (z=-1)
