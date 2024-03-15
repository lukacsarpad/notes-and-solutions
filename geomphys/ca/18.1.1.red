% Maurer-Cartan forms on A(1)

g:= mat((x, y),(0,1));

gi := 1/g;

dg := df(g, x) * dx + df(g, y) * dy;

Om := gi * dg;

OmR := dg * gi;

;;; end;;;