Geometry.CopyMeshingMethod = 1 ;
Mesh.CharacteristicLengthMax = 0.25 ;

NX1  =  21 ; RX1 = 0.95 ; // in front of the cilinder
NX2  =  21 ; RX2 = 1.05 ; // around the cilinder
NX3  =  21 ; RX3 = 1.00 ; // middle of the cilinder
NX4  = 151 ; RX4 = 1.01 ; // behind the cilinder

NY1  =  9 ; RY1  = 1.00 ;
NY1b = 13 ; RY1b = 0.98 ;
NY2  = 31 ; RY2  = 1.05 ;

Z0 =  0.0 ; // plane
SB =  0.1 ; // dummy size
R  =  1.0 ; // unity radius
F  =  2.0 ; // factor around the cilinder

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
Point( 1) = { X0, Y0, Z0, SB};
Point( 2) = { X1, Y0, Z0, SB};
Point( 3) = { X2, Y0, Z0, SB};
Point( 4) = { X5, Y0, Z0, SB};
Point( 5) = { X6, Y0, Z0, SB};
Point( 6) = { X7, Y0, Z0, SB};

Point( 7) = { X3, Y1, Z0, SB};
Point( 8) = { X4, Y1, Z0, SB};

Point( 9) = { X0, Y2, Z0, SB};
Point(10) = {-Y2, Y2, Z0, SB};
Point(11) = { Y2, Y2, Z0, SB};
Point(12) = { X7, Y2, Z0, SB};

Point(13) = { X0, Y3, Z0, SB};
Point(14) = {-Y2, Y3, Z0, SB};
Point(15) = { Y2, Y3, Z0, SB};
Point(16) = { X7, Y3, Z0, SB};

Line(1) = {1,2}; Transfinite Line{1} = NX1 Using Progression RX1 ;
Line(2) = {2,3}; Transfinite Line{2} = NX2 Using Progression 1/RX2 ;
Circle(3)= {3,0,7}; Transfinite Line{3} = NY1 Using Progression RY1 ;
Circle(4)= {7,0,8}; Transfinite Line{4} = NX3 Using Progression RX3 ;
Circle(5)= {8,0,4}; Transfinite Line{5} = NY1b Using Progression 1/RY1b ;
Line(6) = {4,5}; Transfinite Line{6} = NX2 Using Progression RX2 ;
Line(7) = {5,6}; Transfinite Line{7} = NX4 Using Progression RX4 ;

Line( 8) = {9,10}; Transfinite Line{8} = NX1 Using Progression RX1 ;
Circle(9)= {10,0,11};Transfinite Line{9} = NX3 Using Progression RX3 ;
Line(10) = {11,12}; Transfinite Line{10} = NX4 Using Progression RX4 ;

Line(11) = {13,14}; Transfinite Line{11} = NX1 Using Progression RX1 ;
Line(12) = {14,15}; Transfinite Line{12} = NX3 Using Progression RX3 ;
Line(13) = {15,16}; Transfinite Line{13} = NX4 Using Progression RX4 ;

Line(14) = {1,9}; Transfinite Line{14} = NY1 Using Progression RY1 ;
Line(15) = {9,13}; Transfinite Line{15} = NY2 Using Progression RY2 ;
Circle(16)= {2,0,10};Transfinite Line{16} = NY1 Using Progression RY1 ;
Line(17) = {10,14}; Transfinite Line{17} = NY2 Using Progression RY2 ;
Circle(18)= { 5,0,11}; Transfinite Line{18} = NY1b Using Progression RY1b ;
Line(19) = {11,15}; Transfinite Line{19} = NY2 Using Progression RY2 ;
Line(20) = {6,12}; Transfinite Line{20} = NY1b Using Progression RY1b ;
Line(21) = {12,16}; Transfinite Line{21} = NY2 Using Progression RY2 ;
Line(22) = {7,10}; Transfinite Line{22} = NX2 Using Progression RX2 ;
Line(23) = {8,11}; Transfinite Line{23} = NX2 Using Progression RX2 ;

Line Loop(1) = {1,16,-8,-14};
Ruled Surface(1) = {1};
Transfinite Surface(1) = {1,2,10,9};

Line Loop(2) = {8,17,-11,-15};
Ruled Surface(2) = {2};
Transfinite Surface(2) = {9,10,14,13};

Line Loop(3) = {2,3,22,-16};
Ruled Surface(3) = {3};
Transfinite Surface(3) = {2,3,7,10};
Line Loop(4) = {4,23,-9,-22};
Ruled Surface(4) = {4};
Transfinite Surface(4) = {7,8,11,10};
Line Loop(5) = {5,6,18,-23};
Ruled Surface(5) = {5};
Transfinite Surface(5) = {4,5,11,8};

Line Loop(6) = {9,19,-12,-17};
Ruled Surface(6) = {6};
Transfinite Surface(6) = {10,11,15,14};

Line Loop(7) = {7,20,-10,-18};
Ruled Surface(7) = {7};
Transfinite Surface(7) = {5,6,12,11};
Line Loop(8) = {10,21,-13,-19};
Ruled Surface(8) = {8};
Transfinite Surface(8) = {11,12,16,15};

Recombine Surface(1);
Recombine Surface(2);
Recombine Surface(3);
Recombine Surface(4);
Recombine Surface(5);
Recombine Surface(6);
Recombine Surface(7);
Recombine Surface(8);

Symmetry { 0.0,1.0,0.0,0.0 }{Duplicata{Surface{1};}}
Transfinite Line{27} = NX1 Using Progression 1/RX1 ;
Symmetry { 0.0,1.0,0.0,0.0 }{Duplicata{Surface{2};}}
Transfinite Line{32} = NX1 Using Progression 1/RX1 ;
Transfinite Line{33} = NY2 Using Progression 1/RY2 ;
Symmetry { 0.0,1.0,0.0,0.0 }{Duplicata{Surface{3};}}
Symmetry { 0.0,1.0,0.0,0.0 }{Duplicata{Surface{4};}}
Transfinite Line{41} = NX3 Using Progression 1/RX3 ;
Symmetry { 0.0,1.0,0.0,0.0 }{Duplicata{Surface{5};}}
Symmetry { 0.0,1.0,0.0,0.0 }{Duplicata{Surface{6};}}
Transfinite Line{49} = NX3 Using Progression 1/RX3 ;
Symmetry { 0.0,1.0,0.0,0.0 }{Duplicata{Surface{7};}}
Transfinite Line{53} = NX4 Using Progression 1/RX4 ;
Symmetry { 0.0,1.0,0.0,0.0 }{Duplicata{Surface{8};}}
Transfinite Line{57} = NX4 Using Progression 1/RX4 ;


Physical Surface(14) = {1,2,3,4,5,6,7,8,24,29,34,38,42,46,50,54};
Physical Line(201) = {3,4,5,36,39,43};
Physical Line(501) = {14,15,28,33};
Physical Line(801) = {20,21,52,56};
Physical Line(101) = {11,12,13,32,49,57};
