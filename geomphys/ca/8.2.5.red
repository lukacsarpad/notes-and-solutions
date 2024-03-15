z := f(x, y);

rle := {df(z, x) => fx, df(z, y) => fy};
rle2 := {df(z, x, 2) => fxx, df(z, x, y) => fxy, df(z, y, x) => fxy, df(z, y, 2) => fyy};

array x1(3);
x1(1) := 1;
x1(2) := 0;
x1(3) := df(z, x);

array x2(3);
x2(1) := 0;
x2(2) := 1;
x2(3) := df(z, y);

array eps(3, 3, 3);
eps(1, 2, 3) := eps(2, 3, 1) := eps(3, 1, 2) := 1;
eps(3, 2, 1) := eps(2, 1, 3) := eps(1, 3, 2) := -1;

array nv(3);
for i1 := 1:3 do nv(i1) := for i2 := 1:3 sum for i3:=1:3 sum eps(i1, i2, i3) * x1(i2) * x2(i3);

nnv := sqrt(for i1:=1:3 sum nv(i1)^2);

for i1 := 1:3 do nv(i1) := nv(i1) / nnv;

nv(1) where rle;
nv(2) where rle;
nv(3) where rle;

array g(2, 2);

g(1, 1) := for i1 := 1:3 sum x1(i1)*x1(i1);
g(1, 2) := for i1 := 1:3 sum x1(i1)*x2(i1);
g(2, 1) := for i1 := 1:3 sum x2(i1)*x1(i1);
g(2, 2) := for i1 := 1:3 sum x2(i1)*x2(i1);

gdd := mat((g(1, 1), g(1, 2)), (g(2, 1), g(2,2)));

gdd where rle;


guu := 1/gdd;

guu where rle;


array b(2, 2);

b(1, 1) := for i1 := 1:3 sum df(x1(i1), x) * nv(i1);
b(1, 2) := for i1 := 1:3 sum df(x1(i1), y) * nv(i1);
b(2, 1) := for i1 := 1:3 sum df(x2(i1), x) * nv(i1);
b(2, 2) := for i1 := 1:3 sum df(x2(i1), y) * nv(i1);

bdd := mat((b(1,1), b(1, 2)), (b(2, 1), b(2, 2)));

qq := (bdd where rle2)$ (qq where rle); clear qq;

bud := guu * bdd;

qq := (bud where rle2)$ (qq where rle); clear qq;

;;; end ;;;