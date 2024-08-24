load_package "trigsimp";
load_package "linalg";

pauliX := mat((0, 1), (1, 0));

pauliY := mat((0, -i), (i, 0));

pauliZ := mat((1, 0), (0, -1));

hadamard := 1/sqrt(2) * mat((1, 1),(1, -1));

% identity 1

qq1 := hadamard * pauliX * hadamard;

qq2 := pauliZ;

qq1 - qq2;

% identity 2

qq1 := hadamard * pauliY * hadamard;

qq2 := -pauliY;

qq1 - qq2;

% identity 3

qq1 := hadamard * pauliZ * hadamard;

qq2 := pauliX;

qq1 - qq2;

clear qq1, qq2;


;;;end;;;