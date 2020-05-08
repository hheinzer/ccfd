Point(1) = {0, 0, 0, 1.0};
Point(2) = {0.5, 0, 0, 1.0};
Point(3) = {1.5, 0.2309, 0, 1.0};
Point(4) = {1.5, 0.9342, 0, 1.0};
Point(5) = {0, 0.9342, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};

Curve Loop(1) = {5, 1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Curve(301) = {5};
Physical Curve(401) = {3};
Physical Curve(101) = {4, 1, 2};
Physical Surface(1) = {1};

Transfinite Surface {1} = {5, 4, 3, 1};
Transfinite Curve {5, 3} = 201 Using Progression 1;
Transfinite Curve {1} = 51 Using Progression 1;
Transfinite Curve {2} = 101 Using Progression 1;
Transfinite Curve {4} = 151 Using Progression 1;
Recombine Surface {1};
