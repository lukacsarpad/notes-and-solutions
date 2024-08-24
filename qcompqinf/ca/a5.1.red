clear p, q, n, phi, e1, d;

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

procedure rsa(n, e, x); begin scalar retval;
  on modular;
  setmod n;
  retval := x^e;
  off modular;
  return retval;
end;


p := 11;
q := 3;

n := p*q;

phi := (p-1)*(q-1);

% choose e1;

e1 := 7;

gcd(e1, phi);

d := multinv(phi, e1);


% quantum in encoded as
% 18, 22, 2, 15, 21, 22, 14

rsa(n, e1, 18);

rsa(n, e1, 22);

rsa(n, e1, 2);

rsa(n, e1, 15);

rsa(n, e1, 21);

rsa(n, e1, 22);

rsa(n, e1, 14);


% verify

rsa(n, d, rsa(n, e1, 18));
rsa(n, d, rsa(n, e1, 22));
rsa(n, d, rsa(n, e1, 2));
rsa(n, d, rsa(n, e1, 15));
rsa(n, d, rsa(n, e1, 21));
rsa(n, d, rsa(n, e1, 22));
rsa(n, d, rsa(n, e1, 14));

;;;end;;;
