load_package "linalg";

alpha := pi/2^k;
beta := 2*pi/2^k;
gamma := 0;
delta := 0;

pauliX := mat((0, 1), (1, 0));

pauliY := mat((0, -i), (i, 0));

pauliZ := mat((1, 0), (0, -1));

rx := mat((cos(theta/2), -i*sin(theta/2)),(-i*sin(theta/2), cos(theta/2)));

ry := mat((cos(theta/2), -sin(theta/2)), (sin(theta/2), cos(theta/2)));

rz := mat((exp(-i*theta/2), 0), (0, exp(i*theta/2)) );

i2 := make_identity(2);


u := mat((exp(i*(alpha-beta/2-delta/2))*cos(gamma/2), -exp(i*(alpha-beta/2+delta/2))*sin(gamma/2)),
         (exp(i*(alpha+beta/2-delta/2))*sin(gamma/2),  exp(i*(alpha+beta/2+delta/2))*cos(gamma/2)));

rk := mat((1, 0), (0, exp(2*pi*i/2^k)));

rk - u;


A := sub(theta=beta, rz) * sub(theta=gamma/2, ry);

B := sub(theta=-gamma/2, ry) * sub(theta=-(delta+beta)/2, rz);

C := sub(theta=(delta-beta)/2, rz);

A * B * C - i2;

exp(i*alpha) * A * PauliX * B * PauliX * C - rk;





;;;end;;;