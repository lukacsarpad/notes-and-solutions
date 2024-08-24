load_package "linalg";


% constructed from the truth table on page 159
toffoli := mat((1, 0, 0, 0, 0, 0, 0, 0),
               (0, 1, 0, 0, 0, 0, 0, 0),
               (0, 0, 1, 0, 0, 0, 0, 0),
               (0, 0, 0, 1, 0, 0, 0, 0),
               (0, 0, 0, 0, 1, 0, 0, 0),
               (0, 0, 0, 0, 0, 1, 0, 0),
               (0, 0, 0, 0, 0, 0, 0, 1),
               (0, 0, 0, 0, 0, 0, 1, 0));

%invertible?
det(toffoli);

% unitary?
toffoli * hermitian_tp(toffoli) - make_identity(8);

% now check if fig. 4.9 implements it
cnot := mat((1,0,0,0), (0,1,0,0), (0,0,0,1), (0,0,1,0));

hadamard := 1/sqrt(2) * mat((1, 1),(1, -1));

perm12 := mat((1,0,0,0,0,0,0,0),  %000
              (0,0,1,0,0,0,0,0),  %001
              (0,1,0,0,0,0,0,0),  %010
              (0,0,0,1,0,0,0,0),  %011
              (0,0,0,0,1,0,0,0),  %100
              (0,0,0,0,0,0,1,0),  %101
              (0,0,0,0,0,1,0,0),  %110
              (0,0,0,0,0,0,0,1)); %111

perm12 * tp(perm12) - make_identity(8);
perm12 * perm12  - make_identity(8);

perm23 := mat((1,0,0,0,0,0,0,0),  %000
              (0,1,0,0,0,0,0,0),  %001
              (0,0,0,0,1,0,0,0),  %010 -> 100
              (0,0,0,0,0,1,0,0),  %011 -> 101
              (0,0,1,0,0,0,0,0),  %100 -> 010
              (0,0,0,1,0,0,0,0),  %101 -> 011
              (0,0,0,0,0,0,1,0),  %110
              (0,0,0,0,0,0,0,1)); %111

perm23 * tp(perm23) - make_identity(8);
perm23 * perm23 - make_identity(8);

i2 := make_identity(2);

i4 := make_identity(4);

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


tgate := mat((1, 0), (0, exp(i*pi/4)));

sgate := mat((1, 0), (0, i));

t49 := kronecker_product(kronecker_product(tgate, sgate), i2) *
       cnot23 *
       kronecker_product(i2, kronecker_product(hermitian_tp(tgate), i2)) *
       kronecker_product(i4, hadamard) * cnot23 *
       kronecker_product(i2, kronecker_product(hermitian_tp(tgate), tgate)) *
       cnot13 *
       kronecker_product(i4, hermitian_tp(tgate)) *
       cnot12 *
       kronecker_product(i4, tgate) *
       cnot13 *
       kronecker_product(i4, hermitian_tp(tgate)) *
       cnot12 *
       kronecker_product(i4, hadamard);


t49 - toffoli;

;;;end;;;