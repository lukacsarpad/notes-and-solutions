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

S := mat((1,0), (0, i));;

trace(S)/2;

trace(pauliX * S)/2;

trace(pauliY * S)/2;

trace(pauliZ * S)/2;


alpha := pi/4;

theta := -pi/2;

array n(3);
n(1) := 0;
n(2) := 0;
n(3) := 1;

rp := cos(theta/2) * I2 + i * sin(theta/2) * for i1 := 1:3 sum n(i1) * sigma(i1);

trigsimp(S - exp(i*alpha) * rp, trig);


;;; end;;;

