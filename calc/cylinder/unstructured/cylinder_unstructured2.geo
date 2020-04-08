Point(1) = {0., 0., 0., 1.0};
Point(2) = {1., 0., 0., 1.0};
Point(3) = {-1., 0., 0., 1.0};
Point(4) = {11., 0., 0., 1.0};
Point(5) = {-11., 0., 0., 1.0};
Point(6) = {0., 1., 0., 1.0};
Point(7) = {0., 11., 0., 1.0};
Point(8) = {0., -1., 0., 1.0};
Point(9) = {0., -11., 0., 1.0};
Circle(1) = {6, 1, 2};
Circle(2) = {7, 1, 4};
Circle(3) = {2, 1, 8};
Circle(4) = {4, 1, 9};
Circle(5) = {8, 1, 3};
Circle(6) = {9, 1, 5};
Circle(7) = {3, 1, 6};
Circle(8) = {5, 1, 7};
Line Loop(9) = {8, 2, 4, 6};
Plane Surface(10) = {9};
Delete {
  Surface{10};
}
Line Loop(10) = {7, 1, 3, 5};
Plane Surface(11) = {9, 10};
