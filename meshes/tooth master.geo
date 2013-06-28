ms = 0.2;

Point(204) = {-1.3212, 41.2038, 0.0, ms};
Point(205) = {-4.0555, 33.5851, 0.0, ms};

p400 = newp;
Point(p400 + 1) = {-4.059255293888629, 33.99478500466711, 0.0, ms};
Point(p400 + 2) = {-4.028738102321117, 34.40358875365564, 0.0, ms};
Point(p400 + 3) = {-3.97648729198126, 34.81019453455491, 0.0, ms};
Point(p400 + 4) = {-3.905732582075458, 35.21399744684799, 0.0, ms};
Point(p400 + 5) = {-3.819423290338644, 35.61476846791437, 0.0, ms};
Point(p400 + 6) = {-3.719314211619778, 36.01231965675596, 0.0, ms};
Point(p400 + 7) = {-3.606607806019976, 36.40648662857764, 0.0, ms};
Point(p400 + 8) = {-3.482395715070265, 36.79718203771183, 0.0, ms};
Point(p400 + 9) = {-3.347529321878092, 37.18432970274081, 0.0, ms};
Point(p400 + 10) = {-3.202677115821195, 37.567853641935, 0.0, ms};
Point(p400 + 11) = {-3.04842817768492, 37.94769634220152, 0.0, ms};
Point(p400 + 12) = {-2.885051672055972, 38.32370315068622, 0.0, ms};
Point(p400 + 13) = {-2.71296972131127, 38.6958076475121, 0.0, ms};
Point(p400 + 14) = {-2.533266898779081, 39.0642939292798, 0.0, ms};
Point(p400 + 15) = {-2.346421303646304, 39.42920929895138, 0.0, ms};
Point(p400 + 16) = {-2.152200988921003, 39.79025448472555, 0.0, ms};
Point(p400 + 17) = {-1.951586550969515, 40.14778786941439, 0.0, ms};
Point(p400 + 18) = {-1.745735036218295, 40.5023327678182, 0.0, ms};
Point(p400 + 19) = {-1.535476323077682, 40.85428238682979, 0.0, ms};
Spline(400) = {205, p400 + 1, p400 + 2, p400 + 3, p400 + 4, p400 + 5, p400 + 6, p400 + 7, p400 + 8, p400 + 9, p400 + 10, p400 + 11, p400 + 12, p400 + 13, p400 + 14, p400 + 15, p400 + 16, p400 + 17, p400 + 18, p400 + 19, 204};


Point (1) = {-0.9,35,0,4*ms};
Point (2) = {-1.1,34,0,4*ms};

Line(1) = {2, 205};
Line(2) = {204, 1};
Line(3) = {1, 2};

Line Loop(5) = {1,2,3,400};

Plane Surface(1) = {5};
Physical Surface(1) = {1};
Physical Line(1) = {400};   // contact
Physical Line(2) = {3};     // dirichlet
Physical Line(3) = {1,2};     // neumann

