% Examples for problem 11.3.1
load_package "trigsimp";

% sphere of radius r
operator L;

% Laplacian
for all f let L(f) = 1/(a^2 * sin(theta)) * df( sin(theta) * df(f, theta), theta)
  + 1/(a^2 * sin(theta)^2) * df(f, phi, 2);

% coordinates of embedding space
x := a * sin(theta) * cos(phi);
y := a * sin(theta) * sin(phi);
z := a * cos(theta);

% normal vector
nx := sin(theta) * cos(phi);
ny := sin(theta) * sin(phi);
nz := cos(theta);

% sign of the normal vector? 

% the following should all simplify to 0, where H is the mean curvature,
H := - 2 / a;

trigsimp(L(x) - H * nx);

trigsimp(L(y) - H * ny);

trigsimp(L(z) - H * nz);

clear L;
clear x, y, z;
clear nx, ny, nz;


% cylinder
operator L;

for all f let L(f) = 1/a^2 * df(f, phi, 2) + df(f, y, 2);

% coordinates
x := a * sin(phi);
z := a * cos(phi);

% normal vector
nx := sin(phi);
ny := 0;
nz := cos(phi);

% curvature
H := - 1 / a;

trigsimp(L(x) - H * nx);

trigsimp(L(y) - H * ny);

trigsimp(L(z) - H * nz);


clear L;
clear x, y, z;
clear nx, ny, nz;

;;; end ;;;