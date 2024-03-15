I1 := U1 / R;
I2 := U2 / R;

U1 := U1p + Upp;
U2 := U2p + Upp;

array vars2(2), vars1(4);
vars2(1) := I1;
vars2(2) := I2;

vars1(1) := U1p;
vars1(2) := U2p;
vars1(3) := Upp;

vars1(4) := R;

matrix G(2, 4);

for i1 := 1:2 do for i2 := 1:4 do G(i1, i2) := df(vars2(i1), vars1(i2));

matrix V1(4, 4);
V1(1, 1) := sigma1^2;
V1(2, 2) := sigma2^2;
V1(3, 3) := SV^2;
V1(4, 4) := SR^2;

V2 := G * V1 * tp(G);

% covariance
c12 := V2(1, 2);

c12r := coeffn(c12, sr, 2);

c12u := coeffn(c12, sv, 2);

on factor; c12r; off factor;

on factor; c12u; off factor;


c12r - I1 * I2 / R^2;

V2(1, 1) - (sigma1^2 + SV^2 + I1^2 * SR^2)/R^2;
V2(1, 2) - (SV^2 + I1 * I2 * SR^2)/R^2;
V2(2, 2) - (sigma2^2 + SV^2 + I2^2 * SR^2)/R^2;
V2(2, 1) - V2(1, 2);

;;;end;;;
