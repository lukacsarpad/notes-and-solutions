for all r,p,n let prob(r,p,n) = p^r * (1-p)^(n-r)*factorial(n)/factorial(r)/factorial(n-r);

m := 4;
nn := 8;

for all p let qp(p) = (for r := m+1:nn sum prob(r, p, nn));

solpp := num_solve(qp(p) = 0.95, p=0.7);

pp := sub(solpp, p) $

on rounded;
pp;
off rounded;

for all p let qm(p) = (for r := 0:m-1 sum prob(r, p, nn));

solpm := num_solve(qm(p) = 0.95, p=0.2);

pm := sub(solpm, p) $

on rounded;
pm;
off rounded;


;;; end;;;