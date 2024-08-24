load_package "trigsimp";
load_package "linalg";


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

array n1(3), n2(3), n3(3);

n1(1) := n1x;
n1(2) := n1y;
n1(3) := n1z;

n1unit := {n1x^2 + n1y^2 + n1z^2 => 1};

n2(1) := n2x;
n2(2) := n2y;
n2(3) := n2z;

n2unit := {n2x^2 + n2y^2 + n2z^2 => 1};

r1 := rotmat(n1, beta1);

r2 := rotmat(n2, beta2);

s1 := sin(beta1/2);
c1 := cos(beta1/2);

s2 := sin(beta2/2);
c2 := cos(beta2/2);

s3 := sin(beta3/2);
c3 := cos(beta3/2);


for i1 := 1:3 do n3(i1) := (s1 * c2 * n1(i1) + c1 * s2 * n2(i1) + s1 * s2 * for i2 := 1:3 sum for i3 := 1:3 sum eps3(i1, i2, i3) * n2(i2) * n1(i3)) / s3;

r3 := rotmat(n3, beta3);

c12 := c1 * c2 - s1 * s2 * (for i1 := 1:3 sum n1(i1) * n2(i1));

% check norm

qq := for i1 := 1:3 sum n3(i1)^2 $ length(ws);

qq1 := (qq where sin(beta3/2)^2 => (1 - c12^2)) $ length(ws);

qq2 := (qq1 where n1unit) $ length(ws);

qq3 := (qq2 where n2unit) $ length(ws);

qq4 := trigsimp(qq3);

clear qq, qq1, qq2, qq3, qq4;

qq := r2 * r1 - r3 $ length(ws);

qq1 := (qq where {cos(beta3) => cos(beta3/2)^2 - sin(beta3/2)^2, sin(beta3) => 2 * cos(beta3/2) * sin(beta3/2) }) $ length(ws);

qq2 := (qq1 where {cos(beta3/2) => c12, sin(beta3/2)^2 => (1 - c12^2)}) $ length(ws);

qq3 := trigsimp(qq2) $ length(ws);

qq4 := (qq3 where n1unit) $ length(ws);

qq5 := (qq4 where n2unit);

;;; end;;;

