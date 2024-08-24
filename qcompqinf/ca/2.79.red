load_package "linalg";

a2 := mat((1, 1),(1, 1));

e2 := mateigen(a2, mu);

mu21:=sub(solve(first first e2, mu), mu);

mu22:=sub(solve(first second e2, mu), mu);

ev21 := (-third first e2 where arbcomplex(~k) => 1);

n21 := tp(ev21) * ev21;

ev21 := ev21 / sqrt(n21(1, 1));

ev22 := (third second e2 where arbcomplex(~k) => 1);

n22 := tp(ev22) * ev22;

ev22 := ev22 / sqrt(n22(1, 1));

u2 := mat((ev21(1, 1), ev22(1, 1)), (ev21(2, 1), ev22(2, 1)));

d2 := tp(u2) * a2 * u2;

% verify
u2 * d2 * tp(u2) - a2;


%
a3 := mat((1, 1), (1, 0));

% right eigenvectots
e3r := mateigen(a3, mu);

mu31 := sub(first solve(first first e3r, mu), mu);

mu32 := sub(second solve(first first e3r, mu), mu);

d3 := mat((mu31, 0), (0, mu32));

ev3r1 := (sub(mu = mu31, third first e3r) where arbcomplex(~k) => 1);

n3r1 := ev3r1(1, 1)^2 + ev3r1(2, 1)^2;

ev3r1 := ev3r1 / sqrt(n3r1);


ev3r2 := (sub(mu = mu32, third first e3r) where arbcomplex(~k) => 1);

n3r2 := ev3r2(1, 1)^2 + ev3r2(2, 1)^2;

ev3r2 := ev3r2 / sqrt(n3r2);

u3 := mat((ev3r1(1, 1), ev3r2(1, 1)), (ev3r1(2, 1), ev3r2(2, 1)));
% left eigenvectors
e3l := mateigen(tp(a3), mu);

sub(mu = mu31, first first e3l);

sub(mu = mu32, first first e3l);

ev3l1 := (sub(mu = mu31, third first e3r) where arbcomplex(~k) => 1);

n3l1 := ev3l1(1, 1)^2 + ev3l1(2, 1)^2;

ev3l1 := ev3l1 / sqrt(n3r1);


ev3l2 := (sub(mu = mu32, third first e3r) where arbcomplex(~k) => 1);

n3l2 := ev3l2(1, 1)^2 + ev3l2(2, 1)^2;

ev3l2 := ev3l2 / sqrt(n3l2);

v3 := mat((ev3l1(1, 1), ev3l1(2, 1)), (ev3l2(1, 1), ev3l2(2, 1)));

a3 - u3 * d3 * v3;

;;;end;;;