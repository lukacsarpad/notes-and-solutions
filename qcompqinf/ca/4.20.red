load_package "linalg";

had := mat((1, 1), (1, -1))/sqrt(2);
i2 := make_identity(2);

had1 := kronecker_product(i2, had);

had2 := kronecker_product(had, i2);


cz:= mat((1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,-1));

cnot := mat((1,0,0,0), (0,1,0,0), (0,0,0,1), (0,0,1,0));

%     00 01 10 11
% 00  1  0  0  0
% 01  0  0  0  1
% 10  0  0  1  0
% 11  0  1  0  0
cnot2 := mat((1, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 0, 0));


qq := kronecker_product(had, had) * cnot * kronecker_product(had, had);

qq - cnot2;

b0 := mat((1),(0));
b1 := mat((0),(1));

b00 := kronecker_product(b0, b0);
b01 := kronecker_product(b0, b1);
b10 := kronecker_product(b1, b0);
b11 := kronecker_product(b1, b1);

bp := (b0 + b1)/sqrt(2);
bm := (b0 - b1)/sqrt(2);

bpp := kronecker_product(bp, bp);
bpm := kronecker_product(bp, bm);
bmp := kronecker_product(bm, bp);
bmm := kronecker_product(bm, bm);

% some facts
bp - had * b0;

bm - had * b1;

bpp - kronecker_product(had, had) * b00;

bpm - kronecker_product(had, had) * b01;

bmp - kronecker_product(had, had) * b10;

bmm - kronecker_product(had, had) * b11;

% and not the rest of the exercise

kronecker_product(had, had) * cnot * kronecker_product(had, had) * bpp - bpp;

kronecker_product(had, had) * cnot * kronecker_product(had, had) * bpm - bpm;

kronecker_product(had, had) * cnot * kronecker_product(had, had) * bmp - bmm;

kronecker_product(had, had) * cnot * kronecker_product(had, had) * bmm - bmp;


cnot * bpp - bpp;

cnot * bmp - bmp;

cnot * bpm - bmm;

cnot * bmm - bpm;




;;;end;;;