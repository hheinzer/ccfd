SetFactory("OpenCASCADE");
Point(1) = {0, 0, 0, 1.0};
Point(2) = {3.9, 0, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Point(4) = {0.5, 0.9, 0, 1.0};
Point(5) = {0.9, 0.65, 0, 1.0};
Point(6) = {1.1, 0.5, 0, 1.0};
Point(7) = {1.45, 0.4, 0, 1.0};
Point(8) = {1.85, 0.5, 0, 1.0};
Point(9) = {2.05, 0.6, 0, 1.0};
Point(10) = {2.75, 0.85, 0, 1.0};
Point(11) = {3.9, 1, 0, 1.0};

Line(1) = {3, 1};
Line(2) = {1, 2};
Line(3) = {11, 2};
Spline(4) = {3, 4, 5, 6, 7, 8, 9, 10, 11};

Line Loop(1) = {1, 2, -3, -4};
Plane Surface(1) = {1};

Transfinite Line {1, 3} = 50 Using Progression 1.05;
Transfinite Line {4, 2} = 195 Using Progression 1;
Transfinite Surface {1};
Recombine Surface {1};

Physical Surface(1) = {1};
Physical Line(501) = {1};
Physical Line(401) = {3};
Physical Line(101) = {4, 2};
