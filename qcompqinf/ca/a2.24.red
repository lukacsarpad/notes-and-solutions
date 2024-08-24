load_package "linalg";

% S3 group, elements 0...5
% multiplication table
%
%   012345
% 0 012345
% 1 120534
% 2 201453
% 3 345012
% 4 453201
% 5 534120

array multtbl(5,5);
multtbl(0, 0) := 0;
multtbl(0, 1) := 1;
multtbl(0, 2) := 2;
multtbl(0, 3) := 3;
multtbl(0, 4) := 4;
multtbl(0, 5) := 5;

multtbl(1, 0) := 1;
multtbl(1, 1) := 2;
multtbl(1, 2) := 0;
multtbl(1, 3) := 5;
multtbl(1, 4) := 3;
multtbl(1, 5) := 4;

multtbl(2, 0) := 2;
multtbl(2, 1) := 0;
multtbl(2, 2) := 1;
multtbl(2, 3) := 4;
multtbl(2, 4) := 5;
multtbl(2, 5) := 3;

multtbl(3, 0) := 3;
multtbl(3, 1) := 4;
multtbl(3, 2) := 5;
multtbl(3, 3) := 0;
multtbl(3, 4) := 1;
multtbl(3, 5) := 2;

multtbl(4, 0) := 4;
multtbl(4, 1) := 5;
multtbl(4, 2) := 3;
multtbl(4, 3) := 2;
multtbl(4, 4) := 0;
multtbl(4, 5) := 1;

multtbl(5, 0) := 5;
multtbl(5, 1) := 3;
multtbl(5, 2) := 4;
multtbl(5, 3) := 1;
multtbl(5, 4) := 2;
multtbl(5, 5) := 0;

% index of inverses
array invtb(5);
invtb(0) := 0;
invtb(1) := 2;
invtb(2) := 1;
invtb(3) := 3;
invtb(4) := 4;
invtb(5) := 5;



% representation 0, trivial
array rho0(5);

for i1 := 0:5 do rho0(i1) := 1;

% representation 1, sign
array rho1(5);

for i1:= 0:2 do << rho1(i1) := 1; rho1(i1+3) := -1; >>;

% matrix rep
array rho2(5);

rho2(0) := make_identity(2);

rho2(1) := mat((-1, -sqrt(3)), (sqrt(3), -1))/2;

rho2(2) := mat((-1, sqrt(3)), (-sqrt(3), -1))/2;

rho2(3) := mat((-1, 0), (0, 1));

rho2(4) := mat((1, sqrt(3)), (sqrt(3), -1))/2;

rho2(5) := mat((1, -sqrt(3)), (-sqrt(3), -1))/2;

for i1 := 0:5 do for i2 := 0:5 do write rho2(i1) * rho2(i2) - rho2(multtbl(i1, i2));

% characters
array chi0(5), chi1(5), chi2(5);

for i1:=0:5 do <<
    chi0(i1) := rho1(i1);
    chi1(i1) := rho1(i1);
    chi2(i1) := trace(rho2(i1));
>>;

write "Character of the 2d representation:";
for i1:=0:5 do write i1, ", ", chi2(i1);

% orthogonality
for i1 := 0:5 sum chi0(i1)^2;

for i1 := 0:5 sum chi0(i1) * chi1(i1);

for i1 := 0:5 sum chi0(i1) * chi2(i1);

for i1 := 0:5 sum chi1(i1) * chi1(i1);

for i1 := 0:5 sum chi1(i1) * chi2(i1);

for i1 := 0:5 sum chi2(i1) * chi2(i1);


% let us define a function on the group
array f(5);
f(0) := f0;
f(1) := f1;
f(2) := f2;
f(3) := f3;
f(4) := f4;
f(5) := f5;

% calculate the elements of the Fourier transform

% value for rep 0
h0 := 1/sqrt(6) for ig := 0:5 sum f(ig) * rho0(ig);

h1 := 1/sqrt(6) for ig := 0:5 sum f(ig) * rho1(ig);

h2 := sqrt(2)/sqrt(6) for ig := 0:5 sum f(ig) * rho2(ig);

array hatf(5);
hatf(0) := h0;
hatf(1) := h1;
hatf(2) := h2(1, 1);
hatf(3) := h2(1, 2);
hatf(4) := h2(2, 1);
hatf(5) := h2(2, 2);


FT := make_identity(6);
for i1 := 1:6 do for i2:= 1:6 do ft(i1, i2) := coeffn(hatf(i1-1), f(i2-1), 1);

ft;

hermitian_tp(ft)*ft;

% inverse
array fr(5);

for i1:=0:5 do fr(i1) := (hatf(0) * rho0(i1) + hatf(1) * rho1(i1) + sqrt(2) * trace(h2 * rho2(invtb(i1))) ) / sqrt(6);

% check
fr(0);
fr(1);
fr(2);
fr(3);
fr(4);
fr(5);

% read off coefficients
fh2 := mat((fh211, fh212), (fh221, fh222));

for i1:=0:5 do fr(i1) := (fh0 * rho0(i1) + fh1 * rho1(i1) + sqrt(2) * trace(fh2 * rho2(invtb(i1))) ) / sqrt(6);

array fh(5);
fh(0) := fh0;
fh(1) := fh1;
fh(2) := fh211;
fh(3) := fh212;
fh(4) := fh221;
fh(5) := fh222;


IFT := make_identity(6);
for i1 := 1:6 do for i2:= 1:6 do ift(i1, i2) := coeffn(fr(i1-1), fh(i2-1), 1);

ift;

fr(0);
fr(1);
fr(2);
fr(3);
fr(4);
fr(5);

;;;end;;;