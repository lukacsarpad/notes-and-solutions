load_package "trigsimp";
load_package "linalg";


array n(3);
n(1) := cos(alpha) * cos(beta);
n(2) := cos(alpha) * sin(beta);
n(3) := sin(alpha);

array eps3(3,3,3);
eps3(1,2,3) := eps3(2,3,1) := eps3(3,1,2) := 1;
eps3(3,2,1) := eps3(1,3,2) := eps3(2,1,3) := -1;

procedure rotmat(ax, ang); begin matrix rotm(3,3);
  for i1 := 1:3 do for i2 := 1:3 do rotm(i1, i2) := 0;
  for i1 := 1:3 do rotm(i1, i1) := cos(ang);
  for i1 := 1:3 do for i2 := 1:3 do rotm(i1, i2) := rotm(i1, i2) + (1 - cos(ang)) * ax(i1) * ax(i2);
  for i1 := 1:3 do for i2 := 1:3 do rotm(i1, i2) := rotm(i1, i2) + sin(ang) * (for i3 := 1:3 sum eps3(i1, i3, i2) * ax(i3));
  return rotm;
end;

rot3 := rotmat(n, theta);

pauliX := mat((0, 1), (1, 0));
pauliY := mat((0, -i), (i, 0));
pauliZ := mat((1, 0), (0, -1));

array sigma(3);
sigma(1) := pauliX;
sigma(2) := pauliY;
sigma(3) := pauliZ;

I2 := mat((1, 0), (0, 1));

rp := cos(theta/2) * I2 + i * sin(theta/2) * for i1 := 1:3 sum n(i1) * sigma(i1);

rpa := (hermitian_tp(rp) where {repart(alpha)=>alpha, impart(alpha)=>0, repart(beta)=>beta, impart(beta)=>0, repart(theta)=>theta, impart(theta)=>0});

trigsimp(rp*rpa);

array sigmap(3);

for i1 := 1:3 do sigmap(i1) := rp * sigma(i1) * rpa;

array sigmapp(3);

for i1 := 1:3 do sigmapp(i1) := for i2 := 1:3 sum rot3(i1, i2) * sigma(i2);


% are the two equal?
trigsimp(sigmap(1) - sigmapp(1));

trigsimp(sigmap(2) - sigmapp(2));

trigsimp(sigmap(3) - sigmapp(3));

;;; end;;;

