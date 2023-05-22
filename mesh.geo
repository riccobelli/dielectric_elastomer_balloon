dim = 0.005;
td = 0.9;
te = 0.9;
Rm = te;
Ri = Rm*td;
t1 = 1-Rm;
t2 = Rm-Ri;
nthickd = t2*100;
nthicke = t1*100;
ncirc = 101;
Mesh.SaveAll = 1;
//+
Point(1) = {0, -1, 0, 1};
//+
Point(2) = {0, -1+t1, 0, 1};
//+
Point(3) = {0, 0, 0, 1};
//+
Point(4) = {0, 1-t1, 0, 1};
//+
Point(5) = {0, 1, 0, 1};
//+
Point(6) = {1-t1, 0, 0, 1};
//+
Point(7) = {1, 0, 0, 1};
//+
Point(8) = {0, 1-t1-t2, 0, 1};
//+
Point(9) = {0, -1+t1+t2, 0, 1};
//+
Point(10) = {1-t1-t2, 0, 0, 1};
//+
Circle(1) = {1, 3, 7};
//+
Circle(2) = {7, 3, 5};
//+
Circle(3) = {4, 3, 6};
//+
Circle(4) = {6, 3, 2};
Transfinite Line {1,2,3,4} = ncirc;
//+
Line(5) = {5, 4};
//+
Line(6) = {2, 1};
//+
Line(7) = {4, 8};
//+
Line(8) = {2, 9};
Transfinite Line {5,6} = nthicke;
Transfinite Line {7,8} = nthickd;
//+
Circle(9) = {9, 3, 10};
//+
Circle(10) = {10, 3, 8};
Transfinite Line {9,10} = ncirc;
//+
Line(11) = {10, 6};
//+
Line(12) = {6, 7};
Transfinite Line {12} = nthicke;
Transfinite Line {11} = nthickd;
//+
Physical Curve(13) = {7, 10, 11, 3};
//+
Physical Curve(14) = {5, 3, 12, 2};
//+
Physical Curve(15) = {9, 11, 4, 8};
//+
Physical Curve(16) = {6, 1, 12, 4};
//+
Curve Loop(1) = {10, -7, 3, -11};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 3, 12, 2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {4, 6, 1, -12};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {9, 11, 4, 8};
//+
Plane Surface(4) = {4};

Transfinite Surface {1,2,3,4} AlternateLeft;
