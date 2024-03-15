%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% surface of revolution                                                 %
% check metric                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% packages needed
load_package "taylor";
load_package "trigsimp";
load_package "compact";

in "/home/arpi/reduce/r/grfunc.red";
in "/home/arpi/reduce/r/trn.red";

% spatial dimensions
spdim := 2;


% coordinates
operator x;

x(0) := t;
x(1) := x;
x(2) := theta;

%
fp := df(f(x), x);

% metric
array gdd(3,3), guu(3,3);
gdd(0,0) := -1;
gdd(1,1) := (1 + fp^2);
gdd(2,2) := f(x)^2;

operator dxu;
dxu(0) := dt;
dxu(1) := dx;
dxu(2) := dtheta;

ds2 := for j:=0:spdim sum for k:=0:spdim sum gdd(j,k) * dxu(j) * dxu(k);

% calculate inverse metric
comguu();
rdetg:=sqrt(-detg);
rdetg:=sub({abs(-sigm)=sigm,abs(sin(theta))=sin(theta)},rdetg);


%% volume form
filleta();
filleps2();
%filleps4();

comvolform2();

%%%%%%%%%%%%%%%%%%%%%%%%% Christoffel symbols %%%%%%%%%%%%%%%%%%%%%%%%%%%
comchristoffel();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Riemann tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%
comriemann();

%for j:=0:spdim do for k:=0:spdim do for l:=0:spdim do for m:=0:spdim do <<
%    Riedddd(j,k,l,m) := trigsimp(Riedddd(j,k,l,m));
%    Rieuddd(j,k,l,m) := trigsimp(Rieuddd(j,k,l,m));
%    Rieuudd(j,k,l,m) := trigsimp(Rieuudd(j,k,l,m));
%    Rieuuuu(j,k,l,m) := trigsimp(Rieuuuu(j,k,l,m));
%>>;

riemsymm();

bianchialg();

%%%%%%%%%%%%%%%%%%%%%%%% Verify Bianchi identity %%%%%%%%%%%%%%%%%%%%%%%%
bianchidiff();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ricci tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%
comricci();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Weyl tensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comweyl();


gaudd(1,1,1);
gaudd(1,1,2);
gaudd(1,2,2);

gaudd(2,1,1);
gaudd(2,1,2);
gaudd(2,2,2);


% is J = partial / partial theta Killing?
array Ju(2);
Ju(1) := 0;
Ju(2) := 1;

array Jd(2);
for i1 := 1:2 do Jd(i1) := for i2 := 1:2 sum gdd(i1, i2) * Ju(i2);

Jd(1);
Jd(2);

array DJdd(2, 2);

for i1 := 1:2 do for i2 := 1:2 do Djdd(i1, i2) := df(Jd(i1), x(i2)) - for i3 := 1:2 sum gaudd(i3, i2, i1) * jd(i3);

djdd(1, 1);
djdd(1, 2);

djdd(2, 1);
djdd(2, 2);

array kfdd(2, 2);
for i1 := 1:2 do for i2 := 1:2 do kjdd(i1, i2) := djdd(i1, i2) + djdd(i2, i1);

kjdd(1, 1);
kjdd(1, 2);

kjdd(2, 1);
kjdd(2, 2);

% all zero, Ju Killing.

% consider a unit vector, which is to be a tangent to a geodesic.
array exu(2), ethu(2);
exu(1) := 1 / sqrt(gdd(1, 1));
exu(2) := 0;

ethu(1) := 0;
ethu(2) := 1 / sqrt(gdd(2, 2));

for i1 := 1:2 sum for i2 := 1:2 sum gdd(i1, i2) * exu(i1) * exu(i2);

for i1 := 1:2 sum for i2 := 1:2 sum gdd(i1, i2) * ethu(i1) * ethu(i2) where abs(f(x))^2 => f(x)^2;

array tu(2);
for i1 := 1:2 do tu(i1) := sin(alpha) * exu(i1) + cos(alpha) * ethu(i1);

for i1 := 1:2 sum for i2 := 1:2 sum gdd(i1, i2) * tu(i1) * exu(i2);

for i1 := 1:2 sum for i2 := 1:2 sum gdd(i1, i2) * tu(i1) * ethu(i2) where abs(f(x))^2 => f(x)^2;

% angle with J, conserved
for i1 := 1:2 sum for i2 := 1:2 sum gdd(i1, i2) * tu(i1) * Ju(i2) where abs(f(x))^2 => f(x)^2;



end;
