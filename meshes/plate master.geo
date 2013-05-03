/* A square plate */

ms = 0.1;     // Mesh size

Point(1) = {0, -1.01, 0, ms};
Point(2) = {1, -1.01, 0, ms};
Point(3) = {1, -0.01, 0, ms};
Point(4) = {0, -0.01, 0, ms};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = { 5 };
Recombine Surface { 6 };

Physical Surface(1) = {6};       // whole plate
Physical Line(1) = { 1 };        // base
Physical Line(2) = { 2 };        // right side
Physical Line(3) = { 3 };        // top side
Physical Line(4) = { 4 };        // left side

