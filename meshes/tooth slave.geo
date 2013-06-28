ms = 0.2;

Point(9) = {-6.1505, 32.6456, 0.0, ms};
Point(10) = {-2.306, 40.4433, 0.0, ms};

p13 = newp;
Point(p13 + 1) = {-5.894213026551477, 33.00088048089437, 0.0, ms};
Point(p13 + 2) = {-5.641142756322499, 33.35846836838584, 0.0, ms};
Point(p13 + 3) = {-5.39204470841029, 33.71883264488886, 0.0, ms};
Point(p13 + 4) = {-5.147673737821977, 34.08244131839322, 0.0, ms};
Point(p13 + 5) = {-4.908701963763196, 34.449662315622, 0.0, ms};
Point(p13 + 6) = {-4.675643739885761, 34.82067272237829, 0.0, ms};
Point(p13 + 7) = {-4.448993904155743, 35.19562687848379, 0.0, ms};
Point(p13 + 8) = {-4.229157196400831, 35.57461347274177, 0.0, ms};
Point(p13 + 9) = {-4.016407847468043, 35.95762609712398, 0.0, ms};
Point(p13 + 10) = {-3.811009866182843, 36.34465089521157, 0.0, ms};
Point(p13 + 11) = {-3.613290463120258, 36.73568777799213, 0.0, ms};
Point(p13 + 12) = {-3.423862065081771, 37.13079878586794, 0.0, ms};
Point(p13 + 13) = {-3.243417517715011, 37.53006347709396, 0.0, ms};
Point(p13 + 14) = {-3.072639216513933, 37.93355367425733, 0.0, ms};
Point(p13 + 15) = {-2.912068692426598, 38.34123695836122, 0.0, ms};
Point(p13 + 16) = {-2.762146488000203, 38.75300617768742, 0.0, ms};
Point(p13 + 17) = {-2.623622257006458, 39.1688024322027, 0.0, ms};
Point(p13 + 18) = {-2.499160970862679, 39.58887370870776, 0.0, ms};
Point(p13 + 19) = {-2.392155448978624, 40.01358461537234, 0.0, ms};
Spline(13) = {9, p13 + 1, p13 + 2, p13 + 3, p13 + 4, p13 + 5, p13 + 6, p13 + 7, p13 + 8, p13 + 9, p13 + 10, p13 + 11, p13 + 12, p13 + 13, p13 + 14, p13 + 15, p13 + 16, p13 + 17, p13 + 18, p13 + 19, 10};

Point (1) = {-6,40,0,4*ms};

Line(1) = {1, 9};
Line(2) = {10, 1};

Line Loop(5) = {1,2,13};

Plane Surface(1) = {5};
Physical Surface(1) = {1};
Physical Line(1) = {13};  // contact
Physical Line(2) = {2};   // dirichlet
Physical Line(3) = {3};   // neumann
