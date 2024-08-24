A:=mat((1,0), (1,1));

% is it normal?
tp(A)*A - A*tp(A);

j2:=tp(A)*A;

qq := mateigen(j2, mu);

p := first(first(qq));

solp := solve(p=0, mu);

ev := part(first qq, 3);

mu1 := sub(first solp, mu);
ev1 := (sub(first solp, ev) where arbcomplex(~k) => 1);

mu2 := sub(second solp, mu);
ev2 := (sub(second solp, ev) where arbcomplex(~k) => 1);

ev1 := ev1 / sqrt(trace(tp(ev1)*ev1));

ev2 := ev2 / sqrt(trace(tp(ev2)*ev2));

% should be zero
j2 - mu1 * ev1 * tp(ev1) - mu2 * ev2 * tp(ev2);


j := sqrt(mu1) * ev1 * tp(ev1) + sqrt(mu2) * ev2 * tp(ev2);

%
on rounded;
sqrt(5) * j;
off rounded;
%

j := mat((3, 1),(1,2))/sqrt(5);

j2 - j*j;

% A = U J, so U = A J^-1
U := A * (1/J);

A - U * J;

k2 := A * tp(A);

qq := mateigen(k2, mu);

p := first(first(qq));

solp := solve(p=0, mu);

ev := part(first qq, 3);

mu1 := sub(first solp, mu);
ev1 := (sub(first solp, ev) where arbcomplex(~k) => 1);

mu2 := sub(second solp, mu);
ev2 := (sub(second solp, ev) where arbcomplex(~k) => 1);

ev1 := ev1 / sqrt(trace(tp(ev1)*ev1));

ev2 := ev2 / sqrt(trace(tp(ev2)*ev2));

% should be zero
k2 - mu1 * ev1 * tp(ev1) - mu2 * ev2 * tp(ev2);


k := sqrt(mu1) * ev1 * tp(ev1) + sqrt(mu2) * ev2 * tp(ev2);

%
on rounded;
sqrt(5) * k;
off rounded;
%

k := mat((2, 1),(1,3)) / sqrt(5);

k2 - k * k;

A - K*U;

U - (1/K) * A;

;;;end;;;