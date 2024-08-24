load_package("linalg");

had := mat((1, 1), (1, -1))/sqrt(2);

eig := mateigen(had, mu);

la1 := sub(solve(first first eig, mu), mu);

la2 := sub(solve(first second eig, mu), mu);

ev1 := (third first eig where arbcomplex(~x) => 1);

ev2 := (third second eig where arbcomplex(~x) => 1);

had * ev1 - la1 * ev1;

had * ev2 - la2 * ev2;

ev1 * den ev1(1, 1);

ev2 * den ev2(1, 1);


% pedestrian method
p := det(had - mu * make_identity(2));

% same as above
solve(p=0, mu);

p2 := had - la1 * make_identity(2);

m := p2 * den p2(1, 1);

% clearly, a vector orthogonal to the first column is

ev1a := mat((m(1, 2)), (-m(1, 1)));

m * ev1a;

% check if eigenvector with eigenvalue la1
had * ev1a - la1 * ev1a;

p1 := had - la2 * make_identity(2);

m := p1 * den p1(1, 1);

ev2a := mat((-m(1, 2)), (m(1, 1)));

m * ev2a;

% check if eigenvector with eigenvalue la2
had * ev2a - la2 * ev2a;

;;; end;;;