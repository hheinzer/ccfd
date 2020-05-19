//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.6, 0, 0, 1.0};
//+
Point(3) = {0.6, 0.2, 0, 1.0};
//+
Point(4) = {3.0, 0.2, 0, 1.0};
//+
Point(5) = {3.0, 1.0, 0, 1.0};
//+
Point(6) = {0.0, 1.0, 0, 1.0};
//+
Point(7) = {0, 0.2, 0, 1.0};
//+
Point(8) = {0.6, 1, 0, 1.0};
//+
Line(1) = {6, 8};
//+
Line(2) = {8, 5};
//+
Line(3) = {7, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {1, 2};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 1};
//+
Line(8) = {8, 3};
//+
Line(9) = {3, 2};
//+
Line(10) = {5, 4};
//+
Curve Loop(1) = {6, 3, -8, -1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, -9, -3, 7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, -4, -8, 2};
//+
Plane Surface(3) = {3};
//+
Physical Surface(1) = {1, 3, 2};
//+
Physical Curve(301) = {6, 7};
//+
Physical Curve(401) = {10};
//+
Physical Curve(101) = {1, 2, 5};
//+
Physical Curve(201) = {9, 4};
//+
Recombine Surface {2, 1, 3};
//+
Transfinite Surface {2};
//+
Transfinite Surface {1};
//+
Transfinite Surface {3};
//+
Transfinite Curve {10, 8, 6} = 100 Using Progression 1.1;
