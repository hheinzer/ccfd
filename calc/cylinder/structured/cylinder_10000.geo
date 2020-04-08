Geometry.CopyMeshingMethod = 1 ;
Mesh.CharacteristicLengthMax = 0.25 ;

NX1  =  61 ; RX1 = 0.10 ; // in front of the cilinder
NX2  =  41 ; RX2 = 0.10 ; // around the cilinder
NX3  =  61 ; RX3 = 0.10 ; // middle of the cilinder
NX4  = 151 ; RX4 = 0.10 ; // behind the cilinder

NY1  = 30 ; RY1  = 1.00 ;
NY1b = 30 ; RY1b = 0.98 ;
NY2  = 31 ; RY2  = 1.05 ;

Z0 =  0.0 ; // plane
SB =  0.01 ; // dummy size
R  =  1.0 ; // unity radius
F  =  10.0 ; // factor around the cilinder

X0 = -10.0 * R ;
X1 = -(1+F)* R ;
X2 = -1.0 * R ;
X3 = -0.70710678 * R ;
X4 = 0.70710678 * R ;
X5 = 1.0 * R ;
X6 = (1+F)* R ;
X7 = 60.0 * R ;

Y0 = 0.0 ;
Y1 = 0.70710678 * R ;
Y2 = 0.70710678 * R * (1+F) ;
Y3 = 12.0 * R ;

Point( 0) = { 0, 0, Z0, SB};
Point( 2) = { X1, Y0, Z0, SB};
Point( 3) = { X2, Y0, Z0, SB};
Point( 4) = { X5, Y0, Z0, SB};
Point( 5) = { X6, Y0, Z0, SB};

Point( 7) = { X3, Y1, Z0, SB};
Point( 8) = { X4, Y1, Z0, SB};

Point(10) = {-Y2, Y2, Z0, SB};
Point(11) = { Y2, Y2, Z0, SB};

Line(2) = {2,3}; Transfinite Line{2} = NX2 Using Progression 1/RX2 ;
Circle(3)= {3,0,7}; Transfinite Line{3} = NY1 Using Progression RY1 ;
Circle(4)= {7,0,8}; Transfinite Line{4} = NX3 Using Progression RX3 ;
Circle(5)= {8,0,4}; Transfinite Line{5} = NY1b Using Progression 1/RY1b ;
Line(6) = {4,5}; Transfinite Line{6} = NX2 Using Progression RX2 ;

Circle(9)= {10,0,11};Transfinite Line{9} = NX3 Using Progression RX3 ;

Circle(16)= {2,0,10};Transfinite Line{16} = NY1 Using Progression RY1 ;
Circle(18)= { 5,0,11}; Transfinite Line{18} = NY1b Using Progression RY1b ;
Line(22) = {7,10}; Transfinite Line{22} = NX2 Using Progression RX2 ;
Line(23) = {8,11}; Transfinite Line{23} = NX2 Using Progression RX2 ;

Line Loop(3) = {2,3,22,-16};
Ruled Surface(3) = {3};
Transfinite Surface(3) = {2,3,7,10};
Line Loop(4) = {4,23,-9,-22};
Ruled Surface(4) = {4};
Transfinite Surface(4) = {7,8,11,10};
Line Loop(5) = {5,6,18,-23};
Ruled Surface(5) = {5};
Transfinite Surface(5) = {4,5,11,8};

Recombine Surface(3);
Recombine Surface(4);
Recombine Surface(5);

Symmetry { 0.0,1.0,0.0,0.0 }{Duplicata{Surface{3};}}
Symmetry { 0.0,1.0,0.0,0.0 }{Duplicata{Surface{4};}}
Symmetry { 0.0,1.0,0.0,0.0 }{Duplicata{Surface{5};}}


Physical Surface(14) = {3,4,5,24,29,33};
Physical Line(101) = {3,4,5,26,30,34};
Physical Line(301) = {9,16,18,28,32,36};
Transfinite Line {3} = 40 Using Progression 1;
Transfinite Line {4, 3, 26, 30, 34, 5, 16, 9, 18, 36, 32, 28} = 40 Using Progression 1.;
Transfinite Line {27} = 40 Using Progression 1.05;
Transfinite Line {2} = 40 Using Progression 0.95;
Transfinite Line {22} = 40 Using Progression 1.05;
Transfinite Line {23} = 40 Using Progression 1.05;
Transfinite Line {6} = 40 Using Progression 1.05;
Transfinite Line {31} = 40 Using Progression 1.05;
