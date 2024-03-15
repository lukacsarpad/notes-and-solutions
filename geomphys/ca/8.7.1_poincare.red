%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schwarzschild spacetime                                               %
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
x(1) := x1;
x(2) := x2;


% metric
array gdd(3,3), guu(3,3);
gdd(0,0) := -1;
gdd(1,1) := 1/x2^2;
gdd(2,2) := 2/x2^2;

operator dxu;
dxu(0) := dt;
dxu(1) := dth;
dxu(2) := dph;

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


end;
