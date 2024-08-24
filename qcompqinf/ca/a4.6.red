
procedure multinv(n, k); begin scalar i1, retval, qq;
  retval := 0;
  for i1 := 1:n do <<
    on modular;
    setmod n;
    qq := i1 * k;
    if qq = 1 then retval := i1;
    off modular;
  >>;
  return retval;
end;


multinv(24, 17);

;;;end;;;
