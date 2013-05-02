/* A prims extruded from a square */


Point(1) = {0,0,0, 0.1};
Point(2) = {1,0,0, 0.1};
Point(3) = {1,1,0, 0.1};
Point(4) = {0,1,0, 0.1};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

cube[] = Extrude {0,0,-2} { Surface {6}; };


Physical Volume(1)  = {cube[1]}; // cube
Physical Surface(1) = {6};       // base
Physical Surface(2) = {cube[2]}; // side
Physical Surface(3) = {cube[3]}; // side
Physical Surface(4) = {cube[4]}; // side
Physical Surface(5) = {cube[5]}; // side
Physical Surface(6) = {cube[0]}; // bottom (z=-2)
