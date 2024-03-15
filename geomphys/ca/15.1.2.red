procedure vm(x, m);
begin
  array x(4);
  x(1) := m(1, 1);
  x(2) := m(1, 2);
  x(3) := m(2, 1);
  x(4) := m(2, 2);
end;

procedure mv(x);
begin scalar m;
  m:=mat((x(1), x(2)), (x(3), x(4)));
  return m;
end;

% a point in R^4
array p(4);
p(1) := x1;
p(2) := x2;
p(3) := x3;
p(4) := x4;

% corresponding matrix
mp := mv(p);

H := det mp;

array gradH(4);

for i1 := 1:4 do gradH(i1) := df(H, p(i1));

array qq(4);

% test formula
for i1 := 1:4 do qq(i1) := (-1)^( (i1-1)*(4-i1)/2 ) * p(5 - i1);

gradH(1) - qq(1);
gradH(2) - qq(2);
gradH(3) - qq(3);
gradH(4) - qq(4);

clear qq;

% find tangent space
array dp(4);
dp(1) := dx1;
dp(2) := dx2;
dp(3) := dx3;
dp(4) := dx4;

mdp := mv(dp);

tanspace := coeffn(taylortostandard(taylor(det(mp + eps * mdp), eps, 0, 1)), eps, 1);

% check formula
qq := for i1 := 1:4 sum gradH(i1) * dp(i1);

tanspace - qq;
clear qq;

n2gradH := for i1:=1:4 sum gradH(i1)^2;


;;; end;;;