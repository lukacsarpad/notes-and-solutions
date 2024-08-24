load_package "linalg";

cnot12 := mat((1,0,0,0,0,0,0,0),  % 000
              (0,1,0,0,0,0,0,0),  % 001
              (0,0,0,1,0,0,0,0),  % 010 -> 011
              (0,0,1,0,0,0,0,0),  % 011 -> 010
              (0,0,0,0,1,0,0,0),  % 100
              (0,0,0,0,0,1,0,0),  % 101
              (0,0,0,0,0,0,0,1),  % 110 -> 101
              (0,0,0,0,0,0,1,0)); % 111 -> 110

cnot13 := mat((1,0,0,0,0,0,0,0),  % 000
              (0,1,0,0,0,0,0,0),  % 001
              (0,0,1,0,0,0,0,0),  % 010
              (0,0,0,1,0,0,0,0),  % 011
              (0,0,0,0,0,1,0,0),  % 100 -> 101
              (0,0,0,0,1,0,0,0),  % 101 -> 100
              (0,0,0,0,0,0,0,1),  % 110 -> 111
              (0,0,0,0,0,0,1,0)); % 111 -> 110

cnot21 := mat((1,0,0,0,0,0,0,0),  % 000
              (0,0,0,1,0,0,0,0),  % 001 -> 011
              (0,0,1,0,0,0,0,0),  % 010
              (0,1,0,0,0,0,0,0),  % 011 -> 001
              (0,0,0,0,1,0,0,0),  % 100
              (0,0,0,0,0,0,0,1),  % 101 -> 111
              (0,0,0,0,0,0,1,0),  % 110
              (0,0,0,0,0,1,0,0)); % 111 -> 101

cnot23 := mat((1,0,0,0,0,0,0,0),  % 000
              (0,1,0,0,0,0,0,0),  % 001
              (0,0,1,0,0,0,0,0),  % 010
              (0,0,0,1,0,0,0,0),  % 011
              (0,0,0,0,0,0,1,0),  % 100 -> 110
              (0,0,0,0,0,0,0,1),  % 101 -> 111
              (0,0,0,0,1,0,0,0),  % 110 -> 100
              (0,0,0,0,0,1,0,0)); % 111 -> 101

cnot31 := mat((1,0,0,0,0,0,0,0),  % 000
              (0,0,0,0,0,1,0,0),  % 001 -> 101
              (0,0,1,0,0,0,0,0),  % 010
              (0,0,0,0,0,0,0,1),  % 011 -> 111
              (0,0,0,0,1,0,0,0),  % 100
              (0,1,0,0,0,0,0,0),  % 101 -> 001
              (0,0,0,0,0,0,1,0),  % 110
              (0,0,0,1,0,0,0,0)); % 111 -> 011

cnot32 := mat((1,0,0,0,0,0,0,0),  % 000
              (0,1,0,0,0,0,0,0),  % 001
              (0,0,0,0,0,0,1,0),  % 010 -> 110
              (0,0,0,0,0,0,0,1),  % 011 -> 111
              (0,0,0,0,1,0,0,0),  % 100
              (0,0,0,0,0,1,0,0),  % 101
              (0,0,1,0,0,0,0,0),  % 110 -> 010
              (0,0,0,1,0,0,0,0)); % 111 -> 011


toffoli123 := mat((1, 0, 0, 0, 0, 0, 0, 0),  % 000
                  (0, 1, 0, 0, 0, 0, 0, 0),  % 001
                  (0, 0, 1, 0, 0, 0, 0, 0),  % 010
                  (0, 0, 0, 1, 0, 0, 0, 0),  % 011
                  (0, 0, 0, 0, 1, 0, 0, 0),  % 100
                  (0, 0, 0, 0, 0, 1, 0, 0),  % 101
                  (0, 0, 0, 0, 0, 0, 0, 1),  % 110 -> 111
                  (0, 0, 0, 0, 0, 0, 1, 0)); % 111 -> 110

toffoli213 := mat((1,0,0,0,0,0,0,0),  % 000
                  (0,1,0,0,0,0,0,0),  % 001
                  (0,0,1,0,0,0,0,0),  % 010
                  (0,0,0,1,0,0,0,0),  % 011
                  (0,0,0,0,1,0,0,0),  % 100
                  (0,0,0,0,0,0,0,1),  % 101 -> 111
                  (0,0,0,0,0,0,1,0),  % 110
                  (0,0,0,0,0,1,0,0)); % 111 -> 101

toffoli312 := mat((1,0,0,0,0,0,0,0),  % 000
                  (0,1,0,0,0,0,0,0),  % 001
                  (0,0,1,0,0,0,0,0),  % 010
                  (0,0,0,0,0,0,0,1),  % 011 -> 111
                  (0,0,0,0,1,0,0,0),  % 100
                  (0,0,0,0,0,1,0,0),  % 101
                  (0,0,0,0,0,0,1,0),  % 110
                  (0,0,0,1,0,0,0,0)); % 111 -> 011


% partical cyclic permutation
pcp := mat((1,0,0,0,0,0,0,0),  % 000
           (0,0,0,0,0,0,0,1),  % 001 -> 010
           (0,1,0,0,0,0,0,0),  % 010    011
           (0,0,1,0,0,0,0,0),  % 011    100
           (0,0,0,1,0,0,0,0),  % 100    101
           (0,0,0,0,1,0,0,0),  % 101    110
           (0,0,0,0,0,1,0,0),  % 110    111
           (0,0,0,0,0,0,1,0)); % 111 -> 001

% guess logical formulae, we can detect an intput with a
% formula that has xi for xi = 1 and (1 XOR xi) for xi=0, so the
% 8 possible inputs have the following formulae
% 000 (1 Xo x3)(1 Xo x2)(1 Xo x1)  000
% 001 (1 Xo x3)(1 Xo x2) x1        010
% 010 (1 Xo x3) x2      (1 Xo x1)  011
% 011 (1 Xo x3) x2       x1        100
% 100 x3       (1 Xo x2)(1 Xo x1)  101
% 101 x3       (1 Xo x2) x1        110
% 110 x3        x2      (1 Xo x1)  111
% 111 x3        x2       x1        001
% so we read off formulae
% o3 = [(1 Xo x3) x2 x1] + [x3(1 Xo x2)*(1 Xo x1)]
%    + [x3(1 Xo x2)x1] + [x3 x2(1 Xo x1)]
% and (a Xo b)c = ac Xo bc, so
% o3 = [x2 x1 Xo x3 x2 x1] + [x3 Xo x3 x1 Xo x3 x2 Xo x3 x2 x1]
%    + [x3 x2 x1 Xo x3 x1] + [x3 x2 x1 Xo x3 x2]
% which can be simplified:
% note, that the []s are mutually exclusive, can
% XOR'd, and formulae appearing multiple times are cancelled, so
% whar remains is
% o3 = (x2 x1) Xo x3
% note that this is just the output of a Toffoli with x3->o3 being
% the target
% Let us also compute
% o2 = [(1 Xo x3)(1 Xo x2) x1] + [(1 Xo x3) x2      (1 Xo x1)]
%    + [x3       (1 Xo x2) x1] + [x3        x2      (1 Xo x1)]
% which expands to
% o2 = [x1 Xo x2 x1 Xo x3 x1 Xo x3 x2 x1] + [x2 Xo x3 x2 Xo x2 x1 Xo x3 x2 x1]
%    + [x3 x1 Xo x3 x2 x1] + [x3 x2 Xo x3 x2 x1]
% simplified to
% o2 = x1 Xo (x2 x1) Xo (x3 x1) Xo (x3 x2 x1) Xo x2 Xo (x3 x2) Xo (x2 x1) Xo (x3 x2 x1)
%      !     1          2          3             !     4          1          3
%      Xo (x3 x1) Xo (x3 x2 x1) Xo (x3 x2) Xo (x3 x2 x1) = x1 Xo x2,
%         2          3             4          3
% and similarly
% o1 = [(1 Xo x3)x2(1 Xo x1)] + [x3(1 Xo x2)(1 Xo x1)]
%    + [x3 x2(1 Xo x1)] + [x3 x2 x1]
%    = [x2 Xo x3 x2 Xo x2 x1 Xo x3 x2 x1] + [x3 Xo x3 x2 Xo x3 x1 Xo x3 x2 x1]
%       !     1        !        2            !     1        !        2
%    + [x3 x2 Xo x3 x2 x1] + [x3 x2 x1]
%       1        2             2
%    = x2 Xo x2 x1 Xo x3 Xo x3 x1 = x3 Xo x2 Xo x3 x1 Xo x2 x1

% implement with CNOTs and Toffolis
% Note: Toffoli controlling x3 gets o3 right (x3 Xo x1 x2)
%       CNOT                x2 by x1 gets o2 right, remains
% o1 = x3 Xo x2 Xo x3 x1 Xo x2 x1,
%    we can now control by [x3 Xo (x2 x1)] and [x1 Xo x2]
%    we definitely need an x2 and an x3, so we add CNOT for both and see what remains,
%    we get by the CNOTs,
%    x1 Xo (x1 Xo x2) Xo [x3 Xo (x2 x1)] = x2 Xo x3 Xo (x2 x1)
%    1      1      !      !       !
% still missing: x2 Xo x1 x3, which may be added with a Toffoli controlled by (x2 x1) Xo x3 and x1 Xo x2
% product is [(x2 x1) Xo x3][x1 Xo x2] = x3 [x1 Xo x2] and this is XOR'd to the above, yielding
% x2 Xo x3 Xo (x2 x1) Xo x3 [x1 Xo x2] = x2 Xo x3 Xo x2 x1 Xo x3 x1 Xo x3 x2 = x3 Xo x2 Xo x2 x1 Xo x3 x1

u := toffoli123 *
     cnot13 *
     cnot12 *
     cnot21 *
     toffoli312;

u - pcp;



;;;end;;;