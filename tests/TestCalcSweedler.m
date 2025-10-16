%[text] # TestCalcSweedler
% help calcSweedler
[Z,g]=CyclicGroupAlg.getGenerator(3);
ZH=HeisenbergDouble.getGenerator2(Z,'HZ3'); %[output:1a15317c] %[output:704f9238]
S=Z.getSC('antipode');
fa=ZH.make(2,2).split;
gb=ZH.make(3,5).split;

%%
[A,B,C]=calcSweedler('fa{1}|fa{2}') %[output:960c2312] %[output:4ef66d6b] %[output:86e9717a] %[output:6485eb35]

%%
A %[output:11bc8874]
B %[output:443321df]
C %[output:0b105fb7]
A.left %[output:7ee39b67]
A.right %[output:3917b9a6]
B.index %[output:913ea5fc]
%%
%[text] the product on the drinfeld double

str='fa{1}*ev(gb{1}, ev(S^-1,fa{2}{3})*?*fa{2}{1}) | fa{2}{2}*gb{2}';
calcSweedler(str) %[output:2eb1dd8a] %[output:9c6a7655]
'fa{1}*<gb{1}, Sn(fa{2}{3})*?*fa{2}{1},-1)> | fa{2}{2}*gb{2}';
%%
%[text] pentagon equation
'W[1]{1,2}*W[2]{1,3}*W[3]{2,3}=W[4]{2,3}*W[5]{1,2}'
%%
%[text] pentagon equation(coherence)
'<Delta[1]{3},Phi[1]>*<Delta[2]{3},Phi[2]>'
'Phi[1]{2,3,4}*<Delta[1]{2},Phi[2]>*Phi[3]{1,2,3}'
%%
'x{1}+x{2}' %　Δ(x)=x⊗1+1⊗x　+を実装する必要あり
%%
%[text] sweedler notationには複数の意味が混在していて，文脈により解釈が別れるため，解析が非常に難しい．特に，idをパディングするのか，余積をとるのかの記法との判別をつける方法を探したい．例えばR行列など複数同じテンソルを用いる場合は，コピーの何個目かを判別するための添字が必要となる．そのための手段として，\[n\]をつけることでコピー要素の判別を行う．
%[text] R行列は一括して評価が行われるが，テンソルの式の評価の際に，コピーされて評価されるようにする．五角関係式ややんバグスター方程式などではコピーする数分だけ変数を用意するのは面倒である．
%[text] 括弧の種類を区別できたらいいが，通常使うことができる括弧のかずが\[\].{},(),\<\>であり，()は数式の優先順位として用いたい，\<\>はevaluationで用いたい　という要望から，使える括弧が{}と\[\]だけになり，括弧の種類が足りていないということになる．この問題をどうにかできる記法は存在しないだろうか？　文脈で判断させるためには，一回その式の評価を行う必要があり，一貫して同じ評価を行える保証は薄くなる．　これはリスキーである．
%%
%[text] 
A.left.right %[output:2158ccca]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:1a15317c]
%   data: {"dataType":"warning","outputData":{"text":"警告: transpose"}}
%---
%[output:704f9238]
%   data: {"dataType":"warning","outputData":{"text":"警告: transpose"}}
%---
%[output:960c2312]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"3"}}
%---
%[output:4ef66d6b]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"A","value":"     expr: 'fa{1}|fa{2}'\n      num: 1\n     type: 'op'\n     oper: '|'\n     left: [1×1 struct]\n    right: [1×1 struct]"}}
%---
%[output:86e9717a]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"B","value":"    itype: 'dim'\n    index: [2×3 double]\n     rank: 2\n      num: [2 3]\n    value: [1×1 HeisenbergDouble]\n     type: 'vec'"}}
%---
%[output:6485eb35]
%   data: {"dataType":"tabular","outputData":{"columnNames":["type","expr","num","atom","index","kind"],"columns":6,"cornerText":"フィールド","dataTypes":["char","char","double","char","cell","cell"],"header":"フィールドをもつ 1×2 の struct 配列:","name":"C","rows":2,"type":"struct","value":[["'leaf'","'fa{1}'","2","'fa'","1×1 cell","1×1 cell"],["'leaf'","'fa{2}'","3","'fa'","1×1 cell","1×1 cell"]]}}
%---
%[output:11bc8874]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"A","value":"     expr: 'fa{1}|fa{2}'\n      num: 1\n     type: 'op'\n     oper: '|'\n     left: [1×1 struct]\n    right: [1×1 struct]"}}
%---
%[output:443321df]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"B","value":"    itype: 'dim'\n    index: [2×3 double]\n     rank: 2\n      num: [2 3]\n    value: [1×1 HeisenbergDouble]\n     type: 'vec'"}}
%---
%[output:0b105fb7]
%   data: {"dataType":"tabular","outputData":{"columnNames":["type","expr","num","atom","index","kind"],"columns":6,"cornerText":"フィールド","dataTypes":["char","char","double","char","cell","cell"],"header":"フィールドをもつ 1×2 の struct 配列:","name":"C","rows":2,"type":"struct","value":[["'leaf'","'fa{1}'","2","'fa'","1×1 cell","1×1 cell"],["'leaf'","'fa{2}'","3","'fa'","1×1 cell","1×1 cell"]]}}
%---
%[output:7ee39b67]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"ans","value":"     type: 'leaf'\n     expr: 'fa{1}'\n      num: 2\n     atom: 'fa'\n    index: {'1'}\n     kind: {'{'}"}}
%---
%[output:3917b9a6]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"ans","value":"     type: 'leaf'\n     expr: 'fa{2}'\n      num: 3\n     atom: 'fa'\n    index: {'2'}\n     kind: {'{'}"}}
%---
%[output:913ea5fc]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"ans","rows":2,"type":"double","value":[["1","1","1"],["1","2","1"]]}}
%---
%[output:2eb1dd8a]
%   data: {"dataType":"text","outputData":{"text":"nodes = \n  フィールドのない 0×0 の空の <a href=\"matlab:helpPopup('struct')\" style=\"font-weight:bold\">struct<\/a> 配列。\nnodes = \n  フィールドのない 0×0 の空の <a href=\"matlab:helpPopup('struct')\" style=\"font-weight:bold\">struct<\/a> 配列。\n","truncated":false}}
%---
%[output:9c6a7655]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"次を使用中のエラー: <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('assert')\" style=\"font-weight:bold\">assert<\/a>\nrank must be 0 for Nd=1 at atom S^-1\nエラー: <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('calcSweedler>evaluateIndices', '\/Users\/nisimoriyuuya\/Desktop\/programing\/Execution\/algebra\/Core\/tensor\/calcSweedler.m', 192)\" style=\"font-weight:bold\">calcSweedler>evaluateIndices<\/a> (<a href=\"matlab: opentoline('\/Users\/nisimoriyuuya\/Desktop\/programing\/Execution\/algebra\/Core\/tensor\/calcSweedler.m',192,0)\">行 192<\/a>)\n            assert(rank(i)==0, ...\n            ^^^^^^^^^^^^^^^^^^^^^^\nエラー: <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('calcSweedler', '\/Users\/nisimoriyuuya\/Desktop\/programing\/Execution\/algebra\/Core\/tensor\/calcSweedler.m', 120)\" style=\"font-weight:bold\">calcSweedler<\/a> (<a href=\"matlab: opentoline('\/Users\/nisimoriyuuya\/Desktop\/programing\/Execution\/algebra\/Core\/tensor\/calcSweedler.m',120,0)\">行 120<\/a>)\n    [st,nodes]=evaluateIndices(st,nodes,atoms,values);\n    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"}}
%---
%[output:2158ccca]
%   data: {"dataType":"textualVariable","outputData":{"header":"フィールドをもつ struct:","name":"ans","value":"     type: 'leaf'\n    value: 'ev(gb{1},ev(S^-1,fa{2}{3})*?*fa{2}{1})'"}}
%---
