Point(1) = {5,2,0,0.5};

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

Extrude {-10,0,0} {
  Surface{6};
}
