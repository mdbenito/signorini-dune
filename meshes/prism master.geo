/* An extruded prism built with a rectangle inscribed in another one. The mesh is much finer
in the middle*/

ms = 0.1;     // Mesh size
height = 1.0; // Prism height

Point(1) = {-1, -0.01, -1, 4*ms};
Point(2) = {2, -0.01, -1, 4*ms};
Point(3) = {2, -0.01, 2, 4*ms};
Point(4) = {-1, -0.01, 2, 4*ms};

Point(5) = {0, -0.01, 1, ms};
Point(6) = {0, -0.01, 0, ms};
Point(7) = {1, -0.01, 0, ms};
Point(8) = {1, -0.01, 1, ms};

Line(1) = {1, 2};
Line(2) = {3, 2};
Line(3) = {3, 4};
Line(4) = {1, 4};
Line(5) = {6, 1};
Line(6) = {2, 7};
Line(7) = {8, 3};
Line(8) = {4, 5};
Line(9) = {5, 6};
Line(10) = {7, 6};
Line(11) = {7, 8};
Line(12) = {5, 8};

/* Lines with index offset by 100 are the same lines but opposite orientation.
   FIXME: is this actually needed? Shouldn't the physical surfaces be automatically
   aligned properly? */
Line(101) = {2, 1};
Line(102) = {2, 3};
Line(103) = {4, 3};
Line(104) = {4, 1};
Line(105) = {1, 6};
Line(106) = {7, 2};
Line(107) = {3, 8};
Line(108) = {5, 4};
Line(109) = {6, 5};
Line(110) = {6, 7};
Line(111) = {8, 7};
Line(112) = {8, 5};

Line Loop(13) = {1, 6, 10, 5};
Line Loop(14) = {102, 107, 111, 106};
Line Loop(15) = {3, 8, 12, 7};
Line Loop(16) = {104, 105, 109, 108};
Line Loop(17) = {9, 110, 11, 112};

Plane Surface(20) = {13};
Plane Surface(21) = {14};
Plane Surface(22) = {15};
Plane Surface(23) = {16};
Plane Surface(24) = {17};

cube[] = Extrude {0, -height, 0} { 
  Surface {20,21,22,23,24};
  Layers {{2,1,1,1,2},{1/30,1/15,1/10,1/5,1}};
  //Recombine;
};

Physical Volume(1)  = {cube[1], cube[7], cube[13], cube[19], cube[25]}; // whole cube
Physical Surface(1) = {20,21,22,23,24};   // base
Physical Surface(2) = {cube[2]}; // side
Physical Surface(3) = {cube[8]}; // side
Physical Surface(4) = {cube[14]}; // side
Physical Surface(5) = {cube[20]}; // side
Physical Surface(6) = {cube[0], cube[6], cube[12], cube[18], cube[24]}; // bottom (z=-1)
