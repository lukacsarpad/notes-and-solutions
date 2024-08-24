load_package "linalg";
load_package "trigsimp";

% constructed from the truth table on page 159
toffoli := mat((1, 0, 0, 0, 0, 0, 0, 0),
               (0, 1, 0, 0, 0, 0, 0, 0),
               (0, 0, 1, 0, 0, 0, 0, 0),
               (0, 0, 0, 1, 0, 0, 0, 0),
               (0, 0, 0, 0, 1, 0, 0, 0),
               (0, 0, 0, 0, 0, 1, 0, 0),
               (0, 0, 0, 0, 0, 0, 0, 1),
               (0, 0, 0, 0, 0, 0, 1, 0));

cnot13 := mat((1,0,0,0,0,0,0,0),  % 000
              (0,1,0,0,0,0,0,0),  % 001
              (0,0,1,0,0,0,0,0),  % 010
              (0,0,0,1,0,0,0,0),  % 011
              (0,0,0,0,0,1,0,0),  % 100 -> 101
              (0,0,0,0,1,0,0,0),  % 101 -> 100
              (0,0,0,0,0,0,0,1),  % 110 -> 111
              (0,0,0,0,0,0,1,0)); % 111 -> 110

cnot12 := mat((1,0,0,0,0,0,0,0),  % 000
              (0,1,0,0,0,0,0,0),  % 001
              (0,0,0,1,0,0,0,0),  % 010 -> 011
              (0,0,1,0,0,0,0,0),  % 011 -> 010
              (0,0,0,0,1,0,0,0),  % 100
              (0,0,0,0,0,1,0,0),  % 101
              (0,0,0,0,0,0,0,1),  % 110 -> 101
              (0,0,0,0,0,0,1,0)); % 111 -> 110

cnot23 := mat((1,0,0,0,0,0,0,0),  % 000
              (0,1,0,0,0,0,0,0),  % 001
              (0,0,1,0,0,0,0,0),  % 010
              (0,0,0,1,0,0,0,0),  % 011
              (0,0,0,0,0,0,1,0),  % 100 -> 110
              (0,0,0,0,0,0,0,1),  % 101 -> 111
              (0,0,0,0,1,0,0,0),  % 110 -> 100
              (0,0,0,0,0,1,0,0)); % 111 -> 101

rx := mat((cos(theta/2), -i*sin(theta/2)),(-i*sin(theta/2), cos(theta/2)));

ry := mat((cos(theta/2), -sin(theta/2)), (sin(theta/2), cos(theta/2)));

rz := mat((exp(-i*theta/2), 0), (0, exp(i*theta/2)) );

% circuit in ex. 4.26

i2 := make_identity(2);

i4 := make_identity(4);

u := kronecker_product(i4, sub(theta=-pi/4, ry)) *
     cnot12 *
     kronecker_product(i4, sub(theta=-pi/4, ry)) *
     cnot13 *
     kronecker_product(i4, sub(theta=pi/4, ry)) *
     cnot12 *
     kronecker_product(i4, sub(theta=pi/4, ry));

u := trigsimp(u, combine);

% verify on various vectors

b0 := mat((1), (0));

b1 := mat((0), (1));

%
x := mat((x0), (x1));

% value on vectors of the form |0>|0>|x>
v00x := kronecker_product(b0, kronecker_product(b0, x));

u * v00x - toffoli * v00x;

% phase theta(0,0,0) = theta(0,0,1) = 0
theta000 := theta001 := 0;

v01x := kronecker_product(b0, kronecker_product(b1, x));

u * v01x - toffoli * v01x;

% phase theta(0,1,0) = theta(0,1,1) = 0
theta010 := theta011 := 0;

v10x := kronecker_product(b1, kronecker_product(b0, x));

u * v10x;

% phase theta(1,0,0) = 0, theta(1,0,1) = pi
theta100 := 0;
theta101 := pi;

v11x := kronecker_product(b1, kronecker_product(b1, x));

u * v11x - toffoli * v11x;

% phase theta(1,1,0) = theta(1,1,1) = 0;
theta110 := theta111 := 0;

phase := make_identity(8);
phase(1,1) := exp(i*theta000);
phase(2,2) := exp(i*theta001);
phase(3,3) := exp(i*theta010);
phase(4,4) := exp(i*theta011);
phase(5,5) := exp(i*theta100);
phase(6,6) := exp(i*theta101);
phase(7,7) := exp(i*theta110);
phase(8,8) := exp(i*theta111);

phase * u - toffoli;

;;;end;;;
