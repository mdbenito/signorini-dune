/* A square plate */

ms = 0.05;     // Mesh size

Point(1) = {0, 0.01, 0, ms};
Point(2) = {1, 0.01, 0, ms};
Point(3) = {1, 1.01, 0, 4*ms};
Point(4) = {0, 1.01, 0, 4*ms};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = { 5 };
//Recombine Surface { 6 };

Point(10) = {0.2, 0.2, 0, ms};
Point(11) = {0.8, 0.2, 0, ms};
Point{10} In Surface{6};
Point{11} In Surface{6};

Physical Surface(1) = { 6 };     // whole plate
Physical Line(1) = { 1 };        // base
Physical Line(2) = { 2 };        // right side
Physical Line(3) = { 3 };        // top side
Physical Line(4) = { 4 };        // left side

