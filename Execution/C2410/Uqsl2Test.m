% Define
syms q
[O,E,F,K,Ki]=StrUqsl2().getGenerator(q);
A=O.make([3 2],{[] [1]});
B=O.make(3:4,{[] [2 3]});
I=O+1;

%% addition
A+B
C=A-B
A-A
assert(A==A)
O+A
A+3
4+A
%% multiplication
% lb(x,y)=[x,y]はLie bracket
F*E*2
% -K^2*F*E+ K^2*F*E
assert(-K^2*F*E+q/(q^2 - 1)*K -q/(q^2 - 1)*K*K*K +1*E*F*K*K==O)
F^3*K^2*E^2
% カシミール作用素は中心に属する
Omega=F*E+((q - q^-1)^-2)*(q*K-2+q^-1*Ki);
assert(lb(E*F^2*K,Omega)==0)
%% TensorProduct
% テンソル積は縦棒 | で記述する。
% テンソル積 | より+,-の演算子のほうが優先度が高いので注意する。
-(K|F)+(F*E|F)
(F|F)*(E|E)
lb(K|K,E|F)
(1|E)+(E|1)
1|E+E|1 % 1⊗(E+E⊗1)と解釈されてエラーの可能性もある

%% comutiplication

assert(Delta(E*F)-Delta(E)*Delta(F)==0)
assert(Delta(E*K)-Delta(E)*Delta(K)==0)
assert(Delta(F*K)-Delta(F)*Delta(K)==0)
assert(counit(E)==0)
assert(counit(F)==0)

assert(Delta(A*B)-Delta(A)*Delta(B)==0)
assert(counit(E*F+E)==0)
C=B*A;
C2=Delta(C);
C2=C2.lfun_(@fun1);
C2.base=C.base;
C1=C2.calc();
assert(C==C1)

% assert(C1-B*A==0)

function [c,p]=fun1(p)
    c=isempty(find(p{1}<3));
    p=p(2);
end
%% act
v=cellfun(@PolAlg,{1 1},{[1 0] [0 1]});
[v.ctype]=deal("S");
x=v(1);
y=v(2);
x.base.ctype="S";
y.base.ctype="S";

assert(act(E,q*(y|x)-(x|y))==0)
assert(act(F,q*(y|x)-(x|y))==0)
assert(act(F,x|y)==q*(y|y))

%%
