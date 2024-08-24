load_package "trigsimp";
load_package "linalg";

pauliX := mat((0, 1), (1, 0));
pauliY := mat((0, -i), (i, 0));
pauliZ := mat((1, 0), (0, -1));

array sigma(3);
sigma(1) := pauliX;
sigma(2) := pauliY;
sigma(3) := pauliZ;

I2 := mat((1, 0), (0, 1));

Had := mat((1,1),(1,-1))/sqrt(2);

trace(Had);

trace(pauliX * Had);

trace(pauliY * Had);

trace(pauliZ * Had);

array n(3);
n(1) := 1/sqrt(2);
n(2) := 0;
n(3) := 1/sqrt(2);

theta := pi;
alpha := -pi/2;

rp := cos(theta/2) * I2 + i * sin(theta/2) * for i1 := 1:3 sum n(i1) * sigma(i1);

Hadp := exp(i*alpha) * rp;

trigsimp(Had - Hadp);


;;; end;;;

