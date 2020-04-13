Geometry.CopyMeshingMethod = 1 ;
Mesh.CharacteristicLengthMax = 0.25 ;

Z0 = 0.0 ; //plane
NX = 150 ; RX=1.00 ;
NX1 = 25 ;
NX2 = NX-NX1 ;
NY  = 150 ; RY=1.020 ;
SB  =  0.1 ;

Point( 0) = { -1.0, 0.0, Z0, SB};
Point( 1) = {  0.0, 0.0, Z0, SB};
Point( 2) = {  0.0, 5.0, Z0, SB};
Point( 3) = { 5.0, 0.0, Z0, SB};
Point( 4) = { 5.0, 5.0, Z0, SB};
Point( 5) = { -1.0, 5.0, Z0, SB};

Line(1) = {0,1}; Transfinite Line{1} = NX1 Using Progression RX ;
Line(2) = {1,2}; Transfinite Line{2} = NY Using Progression RY ;
Line(3) = {1,3}; Transfinite Line{3} = NX2 Using Progression RX ;
Line(4) = {3,4}; Transfinite Line{4} = NY Using Progression RY ;
Line(5) = {2,4}; Transfinite Line{5} = NX2 Using Progression RX ;
Line(6) = {5,2}; Transfinite Line{6} = NX1 Using Progression RX ;
Line(7) = {0,5}; Transfinite Line{7} = NY Using Progression RY ;

Line Loop(1) = {1,2,-6,-7};
Ruled Surface(1) = {1};
Transfinite Surface {1};
Recombine Surface(1);
Line Loop(2) = {3,4,-5,-2};
Ruled Surface(2) = {2};
Transfinite Surface {2};
Recombine Surface(2);
Physical Surface(14) = {1,2};
Physical Line(101) = {1,6,5};
Physical Line(201) = {3};
Physical Line(301) = {7} ;
Physical Line(401) = {4} ;
