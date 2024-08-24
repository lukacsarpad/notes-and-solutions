% Pauli matrices

array sigma(3);

sigma(0) := mat((1,0),(0,1));
sigma(1) := mat((0,1),(1,0));
sigma(2) := mat((0,-i),(i, 0));
sigma(3) := mat((1,0),(0,-1));

sigma(1) * sigma(2) - sigma(2) * sigma(1) - 2 * i * sigma(3);
sigma(2) * sigma(3) - sigma(3) * sigma(2) - 2 * i * sigma(1);
sigma(3) * sigma(1) - sigma(1) * sigma(3) - 2 * i * sigma(2);

array eps(3, 3, 3);
eps(1, 2, 3) := eps(2, 3, 1) := eps(3, 1, 2) := 1;
eps(3, 2, 1) := eps(2, 1, 3) := eps(1, 3, 2) := -1;

for i1 := 1:3 do for i2 := 1:3 do write sigma(i1) * sigma(i2) - sigma(i2) * sigma(i1) - for i3 := 1:3 sum 2 * i * eps(i1, i2, i3) * sigma(i3);

;;;end;;;
