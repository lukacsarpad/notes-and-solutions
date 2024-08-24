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







;;;end;;;