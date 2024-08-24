load_package "linalg";
load_package "trigsimp";

pauliX := mat((0, 1), (1, 0));
pauliY := mat((0, -i), (i, 0));
pauliZ := mat((1, 0), (0, -1));

hadamard := 1/sqrt(2) * mat((1, 1),(1, -1));

rx := mat((cos(theta/2), -i*sin(theta/2)),(-i*sin(theta/2), cos(theta/2)));
ry := mat((cos(theta/2), -sin(theta/2)), (sin(theta/2), cos(theta/2)));
rz := mat((exp(-i*theta/2), 0), (0, exp(i*theta/2)) );

sub(theta=-pi/2, rx * pauliZ * hermitian_tp(rx)) - pauliY;

sub(theta=-pi/2, rx);

hadamard * pauliZ * hadamard - pauliX;





;;;end;;;