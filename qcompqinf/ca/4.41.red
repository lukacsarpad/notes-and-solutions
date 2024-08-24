load_package "linalg";

psi := mat((psi0), (psi1));

npsi2 :=  hermitian_tp(psi) * psi;

npsi2 := npsi2(1, 1);

npsi := sqrt(npsi2);


rz := mat((exp(-i*theta/2), 0), (0, exp(i*theta/2)) );

i2 := make_identity(2);

% constructed from the truth table on page 159
toffoli := mat((1, 0, 0, 0, 0, 0, 0, 0),
               (0, 1, 0, 0, 0, 0, 0, 0),
               (0, 0, 1, 0, 0, 0, 0, 0),
               (0, 0, 0, 1, 0, 0, 0, 0),
               (0, 0, 0, 0, 1, 0, 0, 0),
               (0, 0, 0, 0, 0, 1, 0, 0),
               (0, 0, 0, 0, 0, 0, 0, 1),
               (0, 0, 0, 0, 0, 0, 1, 0));

hadamard := 1/sqrt(2) * mat((1, 1),(1, -1));

sgate := sgate := mat((1, 0), (0, i));


b0 := mat((1), (0));
b1 := mat((0), (1));


psi3 := kronecker_product(b0, kronecker_product(b0, psi));


ucirc := kronecker_product(hadamard, kronecker_product(hadamard, i2)) *
         toffoli *
         kronecker_product(i2, kronecker_product(i2, sgate)) *
         toffoli *
         kronecker_product(hadamard, kronecker_product(hadamard, i2));

psi32 := ucirc * psi3;

% calculate rz, if cos(theta)=3/5

cth := 3/5;

sth := sqrt(1-cth^2);

% half angles
clear th2, cth2, sth2;

eq1 := 2 * cth2 * sth2 = sth;
eq2 := cth2^2 - sth2^2 = cth;

sol1 := solve(eq1, sth2);

eq2a := sub(sol1, eq2);
eq2b := sub(cth2 = sqrt(x), eq2a);

sol2 := solve(eq2b, x);

th2 := acos(sub(first sol2, sqrt(x)));

cth2 := cos(th2);
sth2:= sin(th2);

clear eq1, eq2, sol1, sol2;

rzth := sub({cos(theta/2) = cth2, sin(theta/2) = sth2}, trigsimp(rz, trig));

psi2 := rzth * psi;

% now, we project ps1 to 0 being qubits 3 and 2, and keep 1

% 000 001 010 011 100 101 110 111
proj := mat((1, 0, 0, 0, 0, 0, 0, 0),
            (0, 1, 0, 0, 0, 0, 0, 0));

psi33 := proj * psi32;

% we note that psi33 has components proportional to psi0, psi1

c0 := coeffn(psi33(1, 1), psi0, 1);
c1 := coeffn(psi33(2, 1), psi1, 1);

psi33 - mat((c0*psi0), (c1*psi1));

% also, it is a phase,

conj(c0) * c0 - conj(c1) * c1;

% probability squared of finding in this state:

prob2 := conj(c0)*c0;

d0 := coeffn((1+i)/sqrt(2) * psi2(1, 1), psi0, 1);

d1 := coeffn((1+i)/sqrt(2) * psi2(2, 1), psi1, 1);

psi2 - mat((d0*psi0), (d1*psi1));

% we need to see of ci = exp(i*alpha) * di

qq0 := c0/d0;

qq1 := c1/d1;


psi33 - sqrt(prob2) * (1+i)/sqrt(2) * psi2;

zgate := mat((1,0), (0, -1));

psi4 := psi32 - kronecker_product(b0, kronecker_product(b0, psi33));

psi4p1 := 1/2/sqrt(2) * i * (1+i)/sqrt(2) * zgate * psi;

psi4p := kronecker_product(-kronecker_product(b0, b1) - kronecker_product(b1, b0) + kronecker_product(b1, b1) , psi4p1);

psi4 - psi4p;


% is the full state vector correct?

psi32 - kronecker_product(kronecker_product(b0, b0), sqrt(prob2) * (1+i)/sqrt(2) * rzth * psi)
    - kronecker_product(-kronecker_product(b0, b1) - kronecker_product(b1, b0) + kronecker_product(b1, b1) , 1/2/sqrt(2) * i * (1+i)/sqrt(2) * zgate * psi);


;;;end;;;