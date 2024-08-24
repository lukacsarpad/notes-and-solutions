load_package "trigsimp";
%load_package "linalg";

array n(3);
n(1) := n1;
n(2) := n2;
n(3) := n3;

nunit := {n1^2 + n2^2 + n3^2 => 1};

operator sigma;
noncom(sigma);

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

rp := cos(theta/2)  + i * sin(theta/2) * for i1 := 1:3 sum n(i1) * sigma(i1);

rpa := sub(theta=-theta, rp);

sigmarules := {sigma(1) * sigma(1) => 1, sigma(1) * sigma(2) => i * sigma(3), sigma(1) * sigma(3) => - i * sigma(2),
               sigma(2) * sigma(1) => -i * sigma(3), sigma(2) * sigma(2) => 1, sigma(2) * sigma(3) => i * sigma(1),
               sigma(3) * sigma(1) => i * sigma(2), sigma(3) * sigma(2) => -i * sigma(1), sigma(3) * sigma(3) => 1};


% is the inverse matrix correct?
trigsimp(((rp * rpa where sigmarules) where nunit), compact);

array pauli(3);
pauli(1) := sigma(1);
pauli(2) := sigma(2);
pauli(3) := sigma(3);

array paulip(3);

for i1 := 1:3 do paulip(i1) := trigsimp((rp * pauli(i1) * rpa where sigmarules), combine);

array paulipp(3);

for i1 := 1:3 do paulipp(i1) := for i2 := 1:3 sum rot3(i1, i2) * pauli(i2);

% do these agree?
paulip(1) - paulipp(1) where nunit;

paulip(2) - paulipp(2) where nunit;

paulip(3) - paulipp(3) where nunit;

;;; end;;;

