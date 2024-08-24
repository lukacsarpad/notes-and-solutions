n := 91;

log2n := logb(n, 2);

on rounded;
log2n;
off rounded;

on rounded;
bmax := floor(log2n);
off rounded;

qq := 0;
for b := 2 : bmax do <<
  x := log2n/b;
  on rounded;
    u1 := floor(2^x);
    u2 := ceiling(2^x);
  off rounded;
  if u1 = u2 then <<
    qq := 1;
    write "N = ", n, " is of the form a^b, with a = ", 2^x, ", b = ", b; 
  >>;
>>;

if qq eq 0 then write "N = ", n, " is not of the form a^b.";


ord := 0;
imax := 10;
x := 4;
for i := 1:imax do begin scalar r;
  on modular;
  setmod n;
  r := x^i;
  write x, "^", i, "= ", r;
  if r = 1 and ord = 0 then ord := i;
  off modular;
end;

write "The order of ", x, " modulo ", n, " is ", ord, ".";

on modular;
  setmod n;
  xr2 := x^(ord/2);
  write "With x = ", x, " and ord = ", ord, ", x^(r/2) = ", xr2, ".";
off modular;

% as this is not -1, we have found a factor

write "Factor found: gcd(x^(r/2) - 1, n) = ", gcd(xr2 - 1, n), ".";

;;;end;;;