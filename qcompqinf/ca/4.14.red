load_package "trigsimp";
load_package "linalg";

pauliX := mat((0, 1), (1, 0));

pauliY := mat((0, -i), (i, 0));

pauliZ := mat((1, 0), (0, -1));

rx := mat((cos(theta/2), -i*sin(theta/2)),(-i*sin(theta/2), cos(theta/2)));

ry := mat((cos(theta/2), -sin(theta/2)), (sin(theta/2), cos(theta/2)));

rz := mat((exp(-i*theta/2), 0), (0, exp(i*theta/2)) );

u := exp(i*alpha) * sub(theta=beta, rz) * sub(theta=gamma, ry) * sub(theta=delta, rz);

up := mat((exp(i*(alpha-beta/2-delta/2))*cos(gamma/2), -exp(i*(alpha-beta/2+delta/2))*sin(gamma/2)),
          (exp(i*(alpha+beta/2-delta/2))*sin(gamma/2),  exp(i*(alpha+beta/2+delta/2))*cos(gamma/2)));

hadamard := 1/sqrt(2) * mat((1, 1),(1, -1));

tgate := mat((1, 0), (0, exp(i*pi/4)));

qq1 := hadamard * tgate * hadamard;

qq2 := sub(theta=pi/4, rx);

qq := trigsimp(qq1/qq2, expon);

qq - mat((exp(i*pi/8), 0), (0, exp(i*pi/8)));

% the relative phase is exp(i*pi/8), i.e.,
% HTH = exp(i*pi/8) * Rz(pi/4);

trigsimp(qq1 - exp(i*pi/8) * qq2, expon);

clear qq1, qq2, qq;

;;;end;;;