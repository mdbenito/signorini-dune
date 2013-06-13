/* A square plate with the top border "adapted"
      +----\-----*-----/-----+
      |     \         /      |
      |      ---------       |
      |                      |
      +----------------------+
*/

ms = 0.1;     // Base mesh size

Point(1) = {-1, -1, 0, 6*ms};
Point(2) = {2, -1, 0, 6*ms};
Point(3) = {2, -0, 0, 3*ms};
Point(4) = {1.5, 0, 0, 2*ms};
Point(5) = {1, -0.3, 0, ms};
Point(6) = {0, -0.3, 0, ms};
Point(7) = {-0.5, 0, 0, 2*ms};
Point(8) = {-1, -0, 0, 3*ms};
Point(9) = {0.5, 0, 0, 0.5*ms};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {7, 9};
Line(10) = {9, 4};

Line Loop(20) = {1, 2, 3, 4, 5, 6, 7, 8};
Line Loop(21) = {4, 5, 6, 9, 10 };

Plane Surface(30) = { 20 };
Plane Surface(31) = { 21 };

Physical Surface(1) = { 30, 31 };    // whole plate
Physical Line(1) = { 1 };            // base
Physical Line(2) = { 2 };            // right side
Physical Line(3) = { 3, 10, 9, 7 };  // top side
Physical Line(4) = { 8 };            // left side

