load_package "linalg";

cnot := mat((1,0,0,0), (0,1,0,0), (0,0,0,1), (0,0,1,0));

rho := mat((rho0000, rho0001, rho0010, rho0011),
           (rho0100, rho0101, rho0110, rho0111),
           (rho1000, rho1001, rho1010, rho1011),
           (rho1100, rho1101, rho1110, rho1111));

cnot*cnot;

rho1 := cnot * rho * cnot;

;;;end;;;