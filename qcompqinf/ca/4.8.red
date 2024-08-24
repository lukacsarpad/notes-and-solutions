load_package "trigsimp";

% any 2 x 2 matrix may be written as

operator sigma;
noncom(sigma);

operator v, w;

U := (alpha + i * beta) + for i1 := 1:3  sum i*(v(i1) + i * w(i1)) * sigma(i1);

relrl := {repart(alpha) => alpha, impart(alpha) => 0, repart(v(~k)) => v(k), impart(v(~k)) => 0, repart(w(~k))=>w(k), impart(w(~k)) =>0};

% Pauli matrices are self adjoint

Ua := (alpha - i * beta) + for i1 := 1:3  sum -i*(v(i1) - i * w(i1)) * sigma(i1);

ua * u;


sigmarules := {sigma(1) * sigma(1) => 1, sigma(1) * sigma(2) => i * sigma(3), sigma(1) * sigma(3) => - i * sigma(2),
               sigma(2) * sigma(1) => -i * sigma(3), sigma(2) * sigma(2) => 1, sigma(2) * sigma(3) => i * sigma(1),
               sigma(3) * sigma(1) => i * sigma(2), sigma(3) * sigma(2) => -i * sigma(1), sigma(3) * sigma(3) => 1};


% u adjoint * u shall be 1
uau := (ua * u where sigmarules);


% coeffn of sigma matrices in Ua * U, all must vanish except 0
array cfn(3);

for i1 := 1:3 do cfn(i1) := coeffn(uau, sigma(i1), 1);

cfn(0) := coeffn(coeffn(coeffn(uau, sigma(1), 0), sigma(2), 0), sigma(3), 0);

uaup := cfn(0) + for i1 := 1:3 sum cfn(i1) * sigma(i1);

uau - uaup;

% the coefficient of the 0th sigma matrix, I2 is
cfn(0);

% this may be written as |alpha|^2 + |beta|^2 + positive terms, so we may write
% alpha + i beta as exp(i*alpha) * cos(theta/2), and the rest must be two vectors, v, w
% such that |v|^2 + |w|^2 = sin(theta/2)^2

array vxw(3);
vxw(1) := v(2) * w(3) - v(3) * w(2);
vxw(2) := v(3) * w(1) - v(1) * w(3);
vxw(3) := v(1) * w(2) - v(2) * w(1);

array cfnp(3);

for i1 := 1:3 do cfnp(i1) := 2 * (beta * v(i1) - alpha * w(i1) - vxw(i1));

cfn(1) - cfnp(1);

cfn(2) - cfnp(2);

cfn(3) - cfnp(3);

% we know that alpha^2 + beta^2 + v^2 + w^2 = 1, so parametrise them
% as alpha = cos(theta/2) * cos(alpha), beta = cos(theta/2) * sin(alpha)
% we shall also notice, that if v, w are not collinear, v, v, v x w 

operator v1;

cfn(1) where {v(~k) => alpha * v1(k), w(~k) => beta * v1(k) };

cfn(2) where {v(~k) => alpha * v1(k), w(~k) => beta * v1(k) };

cfn(3) where {v(~k) => alpha * v1(k), w(~k) => beta * v1(k) };

cfn0 := (cfn(0) where {v(~k) => alpha * v1(k), w(~k) => beta * v1(k) });


;;; end;;;

