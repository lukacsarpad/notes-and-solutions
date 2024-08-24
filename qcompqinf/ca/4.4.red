load_package "trigsimp";
load_package "linalg";

hadamard := mat((1, 1),(1, -1))/sqrt(2);

rx := mat((cos(theta/2), -i*sin(theta/2)),(-i*sin(theta/2), cos(theta/2)));

rz := mat((exp(-i*theta/2), 0), (0, exp(i*theta/2)) );

trigsimp(rx * hermitian_tp(rx) where {impart(theta)=>0, repart(theta)=>theta});

qq:= rx * hadamard;

% can qq be Rz times phase?

trigsimp(qq, expon);

% no

qq := rz * hadamard;

% is this phase times rx?
% no
% conclusion: cannot be two matrices; next: can it be three, e.g. like
% Rx(theta)*Rz(theta')*Rx(theta)? In this case, we rotate by the
% inverse, and get a diagonal matrix

qq := rx * hadamard * rx;

% can this be diagonal?

qq(1, 2);

qq(2, 1);

% yes, if theta/2 = +- pi/4

qq1 := (qq where theta => -pi/2);

i*qq1;

% this looks like rz

trigsimp((-i*qq1 - rz where theta=>pi/2), trig);

% as a result we should get hadamard as

rx1 := (rx where theta => pi/2);
rz1 := (rz where theta => pi/2);

had := i * rx1 * rz1 * rx1;

trigsimp(hadamard - had, trig);

;;;end;;;