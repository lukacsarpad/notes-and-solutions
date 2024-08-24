load_package "linalg";

PauliX := mat((0, 1), (1, 0));
PauliY := mat((0, -i), (i, 0));
PauliZ := mat((1, 0), (0, -1));

i2 := make_identity(2);

array sigma(3);
sigma(0) := i2;
sigma(1) := PauliX;
sigma(2) := PauliY;
sigma(3) := PauliZ;


procedure su2mat(u); begin matrix su2m(3,3);
  begin scalar i1, ua;
    ua := for i1 := 1:3 sum u(i1)^2;
    ua := sqrt(ua);
    su2m := cos(ua/2) * sigma(0) + i * sin(ua/2) / ua * for i1 := 1:3 sum u(i1) * sigma(i1);
  end;
  return su2m;
end;

procedure su2mata(u); begin matrix su2m(3,3);
  begin scalar i1, ua;
    ua := for i1 := 1:3 sum u(i1)^2;
    ua := sqrt(ua);
    su2m := cos(ua/2) * sigma(0) - i * sin(ua/2) / ua * for i1 := 1:3 sum u(i1) * sigma(i1);
  end;
  return su2m;
end;

array x(3), y(3);
x(1) := x1;
x(2) := x2;
x(3) := x3;

y(1) := y1;
y(2) := y2;
y(3) := y3;

xy := for i1 := 1:3 sum x(i1) * y(i1);

xs := for i1 := 1:3 sum x(i1)^2;

ys := for i1 := 1:3 sum y(i1)^2;

xm := sqrt(xs);

ym := sqrt(ys);

ux := su2mat(x);

uy := su2mat(y);

uxa := su2mata(x);

uya := su2mata(y);

m1 := ux - uy;

d := det(m1);

t := trace(m1);


%%% too complicated, other approach: use excer 3.1

m2 := 2 * i2 - uxa * uy - uya * ux;

% this matrix is diagonal, proportional to the identity matrix

m2(1, 2);

m2(2, 1);

m2(1, 1) - m2(2, 2);

m2d := m2(1, 1);

% contains scalar product
coeffn(coeffn(num m2d, x1, 1), y1, 1) - coeffn(coeffn(num m2d, x2, 1), y2, 1);

m2dp := 2 * (1 - cos(xm/2) * cos(ym/2) ) - 2 * sin(xm/2) * sin(ym/2) * xy / xm / ym;;

m2d - m2dp;


m3 := (uxa * uy) + (uya * ux);

m3(1, 2);

m3(2, 1);

m3(1, 1) - m3(2, 2);

m3d := m3(1, 1);

m3dp := 2 * cos(xm/2) * cos(ym/2) + 2 * sin(xm/2) * sin(ym/2) * xy / xm / ym;

m3d - m3dp;

;;;end;;;