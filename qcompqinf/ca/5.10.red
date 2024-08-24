b := 5;
n := 21;

imax := 10;

for i := 1:imax do begin scalar r;
  on modular;
  setmod n;
  r := b^i;
  write b, "^", i, "= ", r;
  off modular;
end;

;;;end;;;