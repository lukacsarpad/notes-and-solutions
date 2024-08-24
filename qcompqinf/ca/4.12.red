load_package "trigsimp";
load_package "linalg";

pauliX := mat((0, 1), (1, 0));

pauliY := mat((0, -i), (i, 0));

pauliZ := mat((1, 0), (0, -1));

rx := mat((cos(theta/2), -i*sin(theta/2)),(-i*sin(theta/2), cos(theta/2)));

ry := mat((cos(theta/2), -sin(theta/2)), (sin(theta/2), cos(theta/2)));

rz := mat((exp(-i*theta/2), 0), (0, exp(i*theta/2)) );

u := exp(i*alpha) * sub(theta=beta, rz) * sub(theta=gamma, ry) * sub(theta=delta, rz);

up := mat((exp(i*(alpha-beta/2-delta/2))*cos(gamma/2), -exp(i*(alpha-beta/2+delta/2))*sin(gamma/2)),
          (exp(i*(alpha+beta/2-delta/2))*sin(gamma/2),  exp(i*(alpha+beta/2+delta/2))*cos(gamma/2)));

u - up;

hadamard := 1/sqrt(2) * mat((1, 1),(1, -1));

hadamardp := sub({alpha=pi/2, beta=0, gamma=pi/2, delta=pi}, u);

trigsimp(hadamard - hadamardp, trig);

A := sub(theta=beta, rz) * sub(theta=gamma/2, ry);

B := sub(theta=-gamma/2, ry) * sub(theta=-(delta+beta)/2, rz);

C := sub(theta=(delta-beta)/2, rz);

A*B*C;

trigsimp(A*B*C);

hadamardpp := exp(i*alpha) * A * pauliX * B * pauliX * C;

hadamardpp := sub({alpha=pi/2, beta=0, gamma=pi/2, delta=pi}, hadamardpp);

trigsimp(trigsimp(hadamard - hadamardpp, expon), trig);

sub({alpha=pi/2, beta=0, gamma=pi/2, delta=pi}, A);

sub({alpha=pi/2, beta=0, gamma=pi/2, delta=pi}, B);

sub({alpha=pi/2, beta=0, gamma=pi/2, delta=pi}, C);

;;;end;;;