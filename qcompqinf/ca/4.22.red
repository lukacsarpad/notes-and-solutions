load_package "linalg";

i2 := make_identity(2);

cnot := mat((1,0,0,0), (0,1,0,0), (0,0,0,1), (0,0,1,0));

ps := mat((1,0), (0, exp(i*alpha)));

psm := mat((1,0), (0, exp(-i*alpha)));

cnot*kronecker_product(psm, i2) * cnot * kronecker_product(ps, ps);

;;;end;;;