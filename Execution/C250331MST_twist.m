g=SweedlerAlg.getGenerator;
HD=HeisenbergDouble.getGenerator2(g,'H(Sweedler)');
[G,W,Wi]=HD.getGW;
% sweedler('G^-2*W{1,1}*G*W{1,2}')
%%
g.verifyHopf(@(x,y)all(eqD(x,y,1e-5),"all"))
%%
Widx=find(~eqD(W.cf,0,1e-5)).';
%%
expr=HD;
for ii=Widx
expr=expr+Wi{ii,2}*G*Wi{ii,1};
end
expr
%%

