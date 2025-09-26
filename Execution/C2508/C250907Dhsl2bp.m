%[text] # Drinfeld Double of Uh(sl2+)
%[text] ## Declare Hopf Algebra
N=5;
[Z,K,E]=Uqsl2BorelSmall.getGenerator(1/N,"L")
Z.verifyHopf(@(x,y)all(eqD(x,y,1e-5),"all"))
%%
E.spec
E*K
K.unit.spec
K.ZERO.spec
Delta(K).spec
K.counit()
E.counit()
z2=E.spec.base.ZERO
z2==Z
%%
s=E.zeros.spec

%%
Usl2Bdual=DualAlg.getGenerator(Z)
%%
N=5;
[Z,K,E]=Uqsl2BorelSmall.getGenerator(1/N,"L",qtype="sym")
Z.verifyHopf(@(x,y)all(isAlways(x==y)))
%%
N=5;
[Z,K,E]=Uqsl2BorelSmall.getGenerator(1/N,"L")
Z.verifyHopf(@(x,y)all(eqD(x,y,1e-5)))

%%
E.sparse
%%

Usl2Bdual.bs










%[text] 
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
