load_package "linalg";

b0 := mat((1), (0));
b1 := mat((0), (1));

pauliX := mat((0, 1), (1, 0));

cnot12 := mat((1,0,0,0),  % 00
              (0,1,0,0),  % 01
              (0,0,0,1),  % 10 -> 11
              (0,0,1,0)); % 11 -> 10

cnot21 := mat((1,0,0,0),  % 00
              (0,0,0,1),  % 01 -> 11
              (0,0,1,0),  % 10
              (0,1,0,0)); % 11 -> 01

i2 := make_identity(2);

rx := mat((cos(theta/2), -i*sin(theta/2)),(-i*sin(theta/2), cos(theta/2)));

rx0 := sub(theta=pi/2, i * rx);

U := kronecker_product(i2, hermitian_tp(rx0)) *
     cnot21 *
     kronecker_product(i2, rx0);


psi1 := u * kronecker_product(b0, b0);

%            00 01 10 11
proj := mat((0, 1, 0, 0),
            (0, 0, 0, 1));

psi2 := proj * psi1;


;;;end;;;