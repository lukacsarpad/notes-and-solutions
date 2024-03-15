U := (1 - 2 * m / r)^(1/2);

array x(3);
x(1) := r;
x(2) := theta;
x(3) := phi;

array gdd(3, 3);

gdd(1, 1) := 1/(1 - 2*m/r);
gdd(2, 2) := r^2;
gdd(3, 3) := r^2 * sin(theta)^2;

h := gdd(1, 1) * gdd(2, 2) * gdd(3, 3);

array gradud(3);

for i1 := 1:3 do gradud(i1) := df(u, x(i1));

qq :=  m/sqrt(1-2*m/r)/r^2;

gradud(1) - qq;
gradud(2);
gradud(3);

clear qq;

array guu(3, 3);

% diag metric
for i1 := 1:3 do guu(i1, i1) := 1 / gdd(i1, i1);

array graduu(3);

for i1 := 1: 3 do graduu(i1) := for i2 := 1:3 sum guu(i1, i2) * gradud(i2);

qq := m /r^2 *sqrt(1-2*m/r);

graduu(1) - qq;
graduu(2);
graduu(3);

clear qq;

array sqhgraduu(3);

for i1 := 1:3 do sqhgraduu(i1) := sqrt(h) * graduu(i1);

shgu1 := (sqhgraduu(1) where sign(sin(theta)) => 1);

nqq := num shgu1;

dqq := den shgu1;

nqq1 := nqq/sin(theta)/m/(r-2*m);

dqq1 := dqq/r;

shgu1a := sqrt(nqq1^2 / dqq1^2) * sin(theta) * m * (r - 2*m) / r;

shgu1b := (shgu1a where {sign(r) => 1, sign(2*m - r) => -1});

lapu := 1/sqrt(h) * for i1 := 1:3 sum df(sqhgraduu(i1), x(i1));


;;; end ;;;