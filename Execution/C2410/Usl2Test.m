% Define
[O,E,F,H]=Usl2.getGenerator;
A=O.make([3 2],{[] [1]});
B=O.make(3:4,{[] [2 3]});
I=O+1;
%% addition
A+B
C=A-B
A-A
O+A
A+3
4+A
%% multiplication
% lb(x,y)=[x,y]はLie bracket
F*E*2
assert(H^2*F-4*F+4*F*H-F*H*H==O)
F^3*H^2*E^2
% カシミール作用素は中心に属する
assert(lb(E^3*F*H*F*E^3,F*E+E*F+H^2*(1/2))==0)

%% TensorProduct
% テンソル積は縦棒 | で記述する。
% テンソル積 | より+,-の演算子のほうが優先度が高いので注意する。
-(H|F)+(F*E|F)
(F|F)*(E|E)
lb(H|H,E|F)
(1|E)+(E|1)
1|E+E|1 % 1⊗(E+E⊗1)と解釈されてエラーの可能性もある

%% comutiplication

assert(Delta(A*B)-Delta(A)*Delta(B)==0)
assert(counit(E*F+E)==0)

C=Delta(B*A);
idx=cellfun(@length,C.pw(:,1))==0;
C2=C.make(C.cf(idx),C.pw(idx,2));
assert(C2-B*A==0)

