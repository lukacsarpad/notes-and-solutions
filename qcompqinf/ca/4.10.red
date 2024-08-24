load_package "trigsimp";
load_package "linalg";

rx := mat((cos(theta/2), -i*sin(theta/2)),(-i*sin(theta/2), cos(theta/2)));

ry := mat((cos(theta/2), -sin(theta/2)), (sin(theta/2), cos(theta/2)));

rz := mat((exp(-i*theta/2), 0), (0, exp(i*theta/2)) );

u := exp(i*alpha) * sub(theta=beta, rz) * sub(theta=gamma, ry) * sub(theta=delta, rz);


;;;end;;;