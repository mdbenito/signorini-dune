/* A plate with a rounded base */

ms = 0.1;     // Mesh size

Point(1) = {0, 0.5, 0, 4*ms};
Point(2) = {0.5, 0, 0, ms};
Point(3) = {1, 0.5, 0, 4*ms};
Point(4) = {1, 1, 0, 6*ms};
Point(5) = {0, 1, 0, 6*ms};
Point(6) = {0.5, 0.5, 0, ms};

Line(1) = {3, 4};
Line(2) = {4, 5};
Line(3) = {5, 1};
Circle(4) = {1, 6, 2};
Circle(5) = {2, 6, 3};

Line Loop(6) = {1, 2, 3, 4, 5};
Plane Surface(7) = { 6 };

Physical Surface(1) = { 7 };     // whole plate
Physical Line(1) = { 4, 5 };     // base
Physical Line(2) = { 1 };        // right side
Physical Line(3) = { 2 };        // top side
Physical Line(4) = { 3 };        // left side


