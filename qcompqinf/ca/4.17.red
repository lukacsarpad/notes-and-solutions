load_package "linalg";

had := mat((1, 1), (1, -1))/sqrt(2);
i2 := make_identity(2);

had1 := kronecker_product(i2, had);

had2 := kronecker_product(had, i2);

u := mat((1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,-1));

cz:= mat((1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,-1));

cnot := mat((1,0,0,0), (0,1,0,0), (0,0,0,1), (0,0,1,0));

%     00 01 10 11
% 00  1  0  0  0
% 01  0  0  0  1
% 10  0  0  1  0
% 11  0  1  0  0
cnot2 := mat((1, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 0, 0));

cnot - had1 * cz * had1;


;;;end;;;