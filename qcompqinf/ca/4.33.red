load_package "linalg";

% computational basis for 1 qubit
b0 := mat((1),(0));
b1 := mat((0),(1));

% computational basis for 2 qubits and their measurement operators
b00 := kronecker_product(b0, b0);
b01 := kronecker_product(b0, b1);
b10 := kronecker_product(b1, b0);
b11 := kronecker_product(b1, b1);

p00 := b00 * hermitian_tp(b00);
p01 := b01 * hermitian_tp(b01);
p10 := b10 * hermitian_tp(b10);
p11 := b11 * hermitian_tp(b11);

% Bell states, eq. (1.23) - (1.26)
bell00 := (b00 + b11)/sqrt(2);
bell01 := (b01 + b10)/sqrt(2);
bell10 := (b00 - b11)/sqrt(2);
bell11 := (b01 - b10)/sqrt(2);

pb00 := bell00 * hermitian_tp(bell00);
pb01 := bell01 * hermitian_tp(bell01);
pb10 := bell10 * hermitian_tp(bell10);
pb11 := bell11 * hermitian_tp(bell11);

% the circuit in ex. 4.33 applies operators before the measurement
% operator, i.e., the measurement is not done on rho, but on U rho
% U^\dagger, so the measurement operators pij are replaced by
% U pij U^dagger

i2 := make_identity(2);

hadamard := hadamard := 1/sqrt(2) * mat((1, 1),(1, -1));

cnot := mat((1,0,0,0), (0,1,0,0), (0,0,0,1), (0,0,1,0));

u := cnot * kronecker_product(hadamard, i2);

ua := hermitian_tp(u);

q00 := u * p00 * ua;
q01 := u * p01 * ua;
q10 := u * p10 * ua;
q11 := u * p11 * ua;

% do these agree? if yes, should all be zero:

q00 - pb00;
q01 - pb01;
q10 - pb10;
q11 - pb11;

;;;end;;;