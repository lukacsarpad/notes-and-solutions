% Pauli matrices

array sigma(3);

sigma(0) := mat((1,0),(0,1));
sigma(1) := mat((0,1),(1,0));
sigma(2) := mat((0,-i),(i, 0));
sigma(3) := mat((1,0),(0,-1));


psi0 := mat((1), (0));
psi1 := mat((0), (1));


psi := exp(i*gamma) * (cos(theta/2) * psi0 + exp(i*phi) * sin(theta/2) * psi1);

psiadj := sub(i=-i, tp(psi));

array v(3);

v0 := trigsimp(psiadj * sigma(0) * psi)$

v(0) := v0(1, 1);


v1 := trigsimp(psiadj * sigma(1) * psi, trig)$

v1 := trigsimp(trigsimp(v1, combine), compact)$

v(1) := v1(1, 1);

v2 := trigsimp(psiadj * sigma(2) * psi, trig)$

v2 := trigsimp(trigsimp(v2, combine), compact)$

v(2) := v2(1, 1);

v3 := trigsimp(psiadj * sigma(3) * psi, trig)$

v3 := trigsimp(trigsimp(v3, combine), compact)$

v(3) := v3(1, 1);

clear v1, v2, v3;


;;;end;;;
