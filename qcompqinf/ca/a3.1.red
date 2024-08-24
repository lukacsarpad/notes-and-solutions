load_package "linalg";

pauliX := mat((0, 1), (1, 0));
pauliY := mat((0, -i), (i, 0));
pauliZ := mat((1, 0), (0, -1));

array sigma(3);
sigma(1) := pauliX;
sigma(2) := pauliY;
sigma(3) := pauliZ;

I2 := make_identity(2);

array n(3);
n(1) := n1;
n(2) := n2;
n(3) := n3;

rn := cos(phi/2) * I2 + i * sin(phi/2) * for i1 := 1:3 sum n(i1) * sigma(i1);

array m(3);
m(1) := m1;
m(2) := m2;
m(3) := m3;

rm := cos(theta/2) * I2 + i * sin(theta/2) * for i1 := 1:3 sum m(i1) * sigma(i1);


mn := for i1 := 1:3 sum m(i1) * n(i1);

unitrl := {n1^2 + n2^2 + n3^2 => 1, m1^2 + m2^2 + m3^2 => 1};

% construct the other norm

x := rm - rn;

xa := sub({theta=-theta, phi=-phi}, x);


x2 := xa * x;

t := trace(x2);

t := (t where unitrl);

tp := 4 * (1 - sin(theta/2)*sin(phi/2)*mn - cos(theta/2)*cos(phi/2));

trigsimp(t-tp);


d := det(x2);

d := trigsimp(d);

d := (d where unitrl);

dp0 := 2 * (1 - sin(theta/2) * sin(phi/2) * mn - cos(theta/2) * cos(phi/2));

dp := dp0^2;

trigsimp(d-dp) where unitrl;



% now consider the matrix x;

t := trace(x);

d := det(x);

d := (d where unitrl);

trigsimp(d-dp0);

disc := t^2/4 - d;

discp := 2 * sin(theta/2) * sin(phi/2) * mn - (sin(theta/2)^2 + sin(phi/2)^2);

trigsimp(disc-discp);


d1 := 2 * sin(theta/2) * sin(phi/2) * s - (sin(theta/2)^2 + sin(phi/2)^2);

discpp := (cos((theta-phi)/2) - cos((theta+phi)/2) ) * mn - 1 + (cos(theta)+cos(phi))/2;

trigsimp(disc - discpp);



;;;end;;;