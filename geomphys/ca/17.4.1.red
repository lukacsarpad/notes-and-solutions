% construct stereographic projection

% upper plane will be z=x+ i y, and according to fig. 1.16, this point is at x, y, 1/2 in R3,

zv := tp(mat((x, y, 1/2)));

% the vector connecting this to the south pole sp is the following

sp := tp(mat((0, 0, -1/2)));

uzv := zv - sp;

% the point on the sphere is the one obtained as

pvp := sp + alpha * uzv;

% and being normed

no := tp(pvp) * pvp;

sol1 := solve(no(1, 1) = 1/4, alpha);

pv := sub(first sol1, pvp);

pp := tp(mat((x, y, (1-x^2-y^2)/2)) / (1 + x^2 + y^2));

pv - pp;

clear pp;

%% similar construction on the lower plane

wv := tp(mat((u, v, -1/2)));

np := tp(mat((0, 0, 1/2)));

uwv := wv - np;

pvp2 := np + alpha * uwv;

no := tp(pvp2) * pvp2;

sol2 := solve(no(1, 1) = 1/4, alpha);

pv2 := sub(first sol2, pvp2);

pp := tp(mat((u, v, (u^2 + v^2 - 1)/2))) / (1 + u^2 + v^2);

pv2 - pp;

clear pp;

% quick verification:

tp(pv)*pv - 1/4;

tp(pv2)*pv2 - 1/4;

% tangent vector, let zeta = xi + i eta and mu = kappa + i nu

tv := df(pv, x) * xi + df(pv, y) * eta;

tv2 := df(pv2, u) * kappa + df(pv2, v) * nu;

% also we know that mu = - zeta / z^2;

tv21 := tp(tv) * tv;

tv21 := tv21(1, 1);

tv22 := tp(tv2) * tv2;

tv22 := tv22(1, 1);

% this yields for the xi^2 + eta^2

zeta22 := v2 / tv21 * (xi^2 + eta^2);

% and for kappa^2 + nu^2

mu2 := v2 / tv22 * (kappa^2 + nu^2);

% invariance if divided by z^2

%
qq1 := zeta22 / (x^2 + y^2);

qq2 := mu2 / (u^2 + v^2);

on factor;
qq1;

qq2;

off factor;

qq1a := sub({x=abs(z), y=0}, qq1);

qq2a := sub({u=1/abs(z), v=0}, qq2);

% do they agree?
qq1a - qq2a;

clear qq1a, qq2a;


% the other expression

qq1 := zeta22 / (1 + x^2 + y^2)^2;

qq2 := mu2 / (1 + u^2 + v^2)^2;

qq1 - qq2;

clear qq1, qq2;


;;; end;;;