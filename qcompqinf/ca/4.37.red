load_package "linalg";

U := 1/2 * mat((1, 1, 1, 1),
               (1, i,-1,-i),
               (1,-1, 1,-1),
               (1,-i,-1, i));

% is it unitary?
hermitian_tp(U) * U - make_identity(4);

% start removing elements

U1 := make_identity(4);

U1(1, 1) := conj(U(1, 1))/sqrt(abs(U(1, 1))^2 + abs(U(2, 1))^2);
U1(1, 2) := conj(U(2, 1))/sqrt(abs(U(1, 1))^2 + abs(U(2, 1))^2);
U1(2, 1) :=      U(2, 1) /sqrt(abs(U(1, 1))^2 + abs(U(2, 1))^2);
U1(2, 2) :=     -U(1, 1) /sqrt(abs(U(1, 1))^2 + abs(U(2, 1))^2);

Up1 := U1 * U;

U2 := make_identity(4);

U2(1, 1) := conj(Up1(1, 1))/sqrt(abs(Up1(1, 1))^2+abs(Up1(3, 1))^2);
U2(1, 3) := conj(Up1(3, 1))/sqrt(abs(Up1(1, 1))^2+abs(Up1(3, 1))^2);
U2(3, 1) :=      Up1(3, 1) /sqrt(abs(Up1(1, 1))^2+abs(Up1(3, 1))^2);
U2(3, 3) :=     -Up1(1, 1) /sqrt(abs(Up1(1, 1))^2+abs(Up1(3, 1))^2);

Up2 := U2 * Up1;

U3 := make_identity(4);

U3(1, 1) := conj(Up2(1, 1))/sqrt(abs(Up2(1, 1))^2+abs(Up2(4, 1))^2);
U3(1, 4) := conj(Up2(4, 1))/sqrt(abs(Up2(1, 1))^2+abs(Up2(4, 1))^2);
U3(4, 1) :=      Up2(4, 1) /sqrt(abs(Up2(1, 1))^2+abs(Up2(4, 1))^2);
U3(4, 4) :=     -Up2(1, 1) /sqrt(abs(Up2(1, 1))^2+abs(Up2(4, 1))^2);

Up3 := U3 * Up2;

% proceed with the next column (2)

U4 := make_identity(4);

d4 := sqrt(abs(Up3(2, 2))^2 + abs(Up3(3, 2))^2);
d4 := (d4 where abs(~z)^2 => z*conj(z));

U4(2, 2) := conj(Up3(2, 2))/d4;
U4(2, 3) := conj(Up3(3, 2))/d4;
U4(3, 2) :=      Up3(3, 2) /d4;
U4(3, 3) :=     -Up3(2, 2) /d4;

Up4 := U4 * Up3;

U5 := make_identity(4);

U5(2, 2) := conj(Up4(2, 2))/sqrt(abs(Up4(2, 2))^2+abs(Up4(4, 2))^2);
U5(2, 4) := conj(Up4(4, 2))/sqrt(abs(Up4(2, 2))^2+abs(Up4(4, 2))^2);
U5(4, 2) :=      Up4(4, 2) /sqrt(abs(Up4(2, 2))^2+abs(Up4(4, 2))^2);
U5(4, 4) :=     -Up4(2, 2) /sqrt(abs(Up4(2, 2))^2+abs(Up4(4, 2))^2);

Up5 := U5 * Up4;

% proceed with third column

U6 := make_identity(4);

U6(3, 3) := conj(Up5(3, 3));
U6(3, 4) := conj(Up5(4, 3));
U6(4, 3) := conj(Up5(3, 4));
U6(4, 4) := conj(Up5(4, 4));

Up6 := U6 * Up5;

V1 := hermitian_tp(U1);
V2 := hermitian_tp(U2);
V3 := hermitian_tp(U3);
V4 := hermitian_tp(U4);
V5 := hermitian_tp(U5);
V6 := hermitian_tp(U6);


% is the decomposition correct?
V1 * V2 * V3 * V4 * V5 * V6 - U;
;;;end;;;