load_package "linalg";

% the two basis vectors on which the operator acts are no. 2 and no. 7, so the Gray code is as follows:
%  010 <- 2
%  011 <- 3
%  111 <- 7

% we construct the following C^2 NOT operator

cnot12n3 := mat((1,0,0,0,0,0,0,0),  % 000
                (0,1,0,0,0,0,0,0),  % 001
                (0,0,0,1,0,0,0,0),  % 010
                (0,0,1,0,0,0,0,0),  % 011
                (0,0,0,0,1,0,0,0),  % 100
                (0,0,0,0,0,1,0,0),  % 101
                (0,0,0,0,0,0,1,0),  % 110
                (0,0,0,0,0,0,0,1)); % 111

% the following controlled-U acts on qubit 3 if qubits 2 and 1 are 1 1

utilde := mat((1,0,0,0,0,0,0,0),  % 000
              (0,1,0,0,0,0,0,0),  % 001
              (0,0,1,0,0,0,0,0),  % 010
              (0,0,0,a,0,0,0,c),  % 011
              (0,0,0,0,1,0,0,0),  % 100
              (0,0,0,0,0,1,0,0),  % 101
              (0,0,0,0,0,0,1,0),  % 110
              (0,0,0,b,0,0,0,d)); % 111

% the gate to be implemented is

u := mat((1,0,0,0,0,0,0,0),
         (0,1,0,0,0,0,0,0),
         (0,0,a,0,0,0,0,c),
         (0,0,0,1,0,0,0,0),
         (0,0,0,0,1,0,0,0),
         (0,0,0,0,0,1,0,0),
         (0,0,0,0,0,0,1,0),
         (0,0,b,0,0,0,0,d));


% verify unitarity
cnot12n3 * hermitian_tp(cnot12n3) - make_identity(8);

% verify hermiticity
cnot12n3 - hermitian_tp(cnot12n3);

% verify implementation

cnot12n3 * utilde * cnot12n3 - u;


;;;end;;;