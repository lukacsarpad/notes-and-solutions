load_package "linalg";

i2 := make_identity(2);

% elementary gates
pauliX := mat((0, 1), (1, 0));

pauliY := mat((0, -i), (i, 0));

pauliZ := mat((1, 0), (0, -1));

rx := mat((cos(theta/2), -i*sin(theta/2)),(-i*sin(theta/2), cos(theta/2)));

ry := mat((cos(theta/2), -sin(theta/2)), (sin(theta/2), cos(theta/2)));

rz := mat((exp(-i*theta/2), 0), (0, exp(i*theta/2)) );

% x1 is controlled by 2
cnot12 := mat((1,0,0,0), (0,1,0,0), (0,0,0,1), (0,0,1,0));

% x2 is controlled by 1
cnot21 := mat((1,0,0,0),  % 00
              (0,0,0,1),  % 01
              (0,0,1,0),  % 10
              (0,1,0,0)); % 11



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 := kronecker_product(i2, paulix);
x2 := kronecker_product(paulix, i2);

cnot21 * x1 * cnot21 - x1 * x2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y1 := kronecker_product(i2, pauliy);

cnot21 * y1 * cnot21 - y1 * x2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z1 := kronecker_product(i2, pauliz);

cnot21 * z1 * cnot21 - z1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnot21 * x2 * cnot21 - x2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y2 := kronecker_product(pauliy, i2);

cnot21 * y2 * cnot21 - z1 * y2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z2 := kronecker_product(pauliz, i2);

cnot21 * z2 * cnot21 - z1 * z2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rz1 := kronecker_product(i2, rz);

rz1 * cnot21 - cnot21 * rz1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx2 := kronecker_product(rx, i2);

rx2 * cnot21 - cnot21 * rx2;

;;;end;