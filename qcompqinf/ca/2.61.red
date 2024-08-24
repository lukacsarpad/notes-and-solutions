load_package "linalg";

array sigma(3);

sigma(1) := mat((0, -i), (i, 0));
sigma(2) := mat((0, 1), (1, 0));
sigma(3) := mat((1, 0), (0, -1));

array vecv(3);
v(1) := v1;
v(2) := v2;
v(3) := v3;

vsigma := for i1 := 1:3 sum v(i1) * sigma(i1);

eig := mateigen(vsigma, mu);

eig1 := (eig where v1^2 + v2^2 + v3^2 => 1);

w1 := (third first eig1 where {mu =>1, arbcomplex(~x)=>1});

w2 := (third first eig1 where {mu =>-1, arbcomplex(~x)=>1});

vsigma * w1 - w1 where v1^2 + v2^2 + v3^2 => 1;

vsigma * w2 + w2 where v1^2 + v2^2 + v3^2 => 1;

b1 := mat((1), (0));

b2 := mat((0), (1));

p1 := (make_identity(2) + vsigma)/2;

p2 := (make_identity(2) - vsigma)/2;

p1^2 - p1 where v1^2 + v2^2 + v3^2 => 1;

p2^2 - p2 where v1^2 + v2^2 + v3^2 => 1;


p1 * w1 - w1 where v1^2 + v2^2 + v3^2 => 1;

p1 * w2 where v1^2 + v2^2 + v3^2 => 1;

p2 * w1 where v1^2 + v2^2 + v3^2 => 1;

p2 * w2 - w2 where v1^2 + v2^2 + v3^2 => 1;

psi := p1 * b1;

psin := hermitian_tp(psi) * psi;

psin := (psin where {repart(v1) => v1, impart(v1) => 0});
psin := (psin where {repart(v2) => v2, impart(v2) => 0});
psin := (psin where {repart(v3) => v3, impart(v3) => 0});
psin := (psin where v1^2 + v2^2 + v3^2 => 1);

psi := psi / sqrt(psin(1, 1));





;;; end;;;