q1 := int(1/y^2, y,x,x+1);

q2 := 2/3/x^2;

qq := q1-q2;

qqp := (x/3-2/3)/x^2/(x+1);

qq - qqp;

% so qq is increasing, and positive if x is greater than
solve(num qq = 0, x);

% verification:
sign(qq) where {sign(x-2)=>1, sign(x) =>1};

;;;end;;;