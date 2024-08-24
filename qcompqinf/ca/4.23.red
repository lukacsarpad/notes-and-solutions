load_package "trigsimp";
load_package "linalg";

pauliX := mat((0, 1), (1, 0));

pauliY := mat((0, -i), (i, 0));

pauliZ := mat((1, 0), (0, -1));

rx := mat((cos(theta/2), -i*sin(theta/2)),(-i*sin(theta/2), cos(theta/2)));

ry := mat((cos(theta/2), -sin(theta/2)), (sin(theta/2), cos(theta/2)));

rz := mat((exp(-i*theta/2), 0), (0, exp(i*theta/2)) );

i2 := make_identity(2);

% parameters for paulix

u := mat((exp(i*(alpha-beta/2-delta/2))*cos(gamma/2), -exp(i*(alpha-beta/2+delta/2))*sin(gamma/2)),
         (exp(i*(alpha+beta/2-delta/2))*sin(gamma/2),  exp(i*(alpha+beta/2+delta/2))*cos(gamma/2)));

sub({gamma=pi, alpha=pi/2, beta=0, delta=pi}, u) - paulix;

% parameters for rx

sub({gamma=theta, alpha=0, beta=-pi/2, delta=pi/2}, u) -rx;


% decomposition of rx

A0 := sub(theta=beta, rz) * sub(theta=gamma/2, ry);

B0 := sub(theta=-gamma/2, ry) * sub(theta=-(delta+beta)/2, rz);

C0 := sub(theta=(delta-beta)/2, rz);

A := sub({gamma=theta, alpha=0, beta=-pi/2, delta=pi/2}, A0);
B := sub({gamma=theta, alpha=0, beta=-pi/2, delta=pi/2}, B0);
C := sub({gamma=theta, alpha=0, beta=-pi/2, delta=pi/2}, C0);

trigsimp(A0*B0*C0);

trigsimp(A*pauliX*B*pauliX*C - rx);

cnot := mat((1,0,0,0), (0,1,0,0), (0,0,0,1), (0,0,1,0));

% now construct the full controlled rz
crx := kronecker_product(i2, a) * cnot * kronecker_product(i2, b) * cnot * kronecker_product(i2, c);

crx := trigsimp(crx, combine);

% can we simplify? does something commute?

trigsimp(kronecker_product(i2, a) * cnot - cnot * kronecker_product(i2, a), trig);

trigsimp(kronecker_product(i2, b) * cnot - cnot * kronecker_product(i2, b), trig);

trigsimp(kronecker_product(i2, c) * cnot - cnot * kronecker_product(i2, c), trig);

%
D := 1/C;

trigsimp(A * pauliX * B * pauliX * C - Rx, trig);

trigsimp(A * pauliX * B * (PauliX * C * PauliX) * PauliX - Rx, trig);

trigsimp(A * pauliX * C * D * B * (PauliX * C * PauliX) * PauliX - Rx, trig);

trigsimp(A * (PauliX * C * PauliX) * pauliX * D * B * (PauliX * C * PauliX) * PauliX - Rx, trig);

A1 := A * (PauliX * C * PauliX);

B1 := D * B * (PauliX * C * PauliX);

qq := kronecker_product(i2, a1) * cnot * kronecker_product(i2, b1) * cnot;
% * kronecker_product(i2, c);


;;;end;;;