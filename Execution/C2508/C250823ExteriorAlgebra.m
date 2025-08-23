%[text] # Exterior Algebra Example
%[text] First, define StrExtAlg Class, which is derived from StrAlg, 
%[text] Exterior algebra is not a Hopf algebra, but a Hopf monoid with super vector space.
%[text] in a multiplication of two higher rank tensors, Z\_2 graded is needed.
N=3;
[Z,X]=StrExtAlg.getGenerator(N) %[output:7d629f1d] %[output:70f1d02e]
eval(join(arrayfun(@(i)sprintf("X%d=X(%d);",i,i),1:N)));
%%
X1*X2 %[output:2d29d69a]
X2*X1 %[output:68e9b3f7]
X1^2 %[output:84490f4d]
X1*X1 %[output:9a1a7a84]
X3*X1*X2 %[output:9ac4dab0]
(X1|X2)*(X2|X3) %[output:6464cfa9]
(X3|X2)*(X2|X3) %[output:3a9438be]
(X1|1)*(1|X1) %[output:5b67a257]
(1|X1)*(X1|1) %[output:0d57966f]
Delta(X1*X2) %[output:0fca7793]
Delta(X2*X2) %[output:12d298ff]
Delta(X1*X2*X3) %[output:6c0077b2]
%%
[Z,X]=VectExtAlg.getGenerator(3) %[output:0f8e0a60] %[output:1ed56b43]
eval(join(arrayfun(@(i)sprintf("X%d=X(%d);",i,i),1:N)));
X1*X2*X1 %[output:53bfa2a8]
X1*X2 %[output:81504f67]
X2*X3*X1 %[output:01c915eb]
X2*X1*X3 %[output:5a3c6adc]
%%
Delta(X1*X2) %[output:9222b755]
Delta(X1*X2*X3) %[output:8083797b]



%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:7d629f1d]
%   data: {"dataType":"textualVariable","outputData":{"name":"Z","value":"(0)*1"}}
%---
%[output:70f1d02e]
%   data: {"dataType":"text","outputData":{"text":"X =\n  3 個のエントリをもつ <a href=\"matlab:helpPopup('dictionary')\" style=\"font-weight:bold\">dictionary<\/a> (<strong>double<\/strong> ⟼ <strong>StrAlg<\/strong>):\n    1 ⟼ 1×1 StrExtAlg\n    2 ⟼ 1×1 StrExtAlg\n    3 ⟼ 1×1 StrExtAlg\n","truncated":false}}
%---
%[output:2d29d69a]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"(1)*X1*X2"}}
%---
%[output:68e9b3f7]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"(-1)*X1*X2"}}
%---
%[output:84490f4d]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"(1)*X1*X1"}}
%---
%[output:9a1a7a84]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"(0)*1"}}
%---
%[output:9ac4dab0]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"(1)*X3*X1*X2"}}
%---
%[output:6464cfa9]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"(-1)*(X1*X2|X2*X3)"}}
%---
%[output:3a9438be]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"(1)*(X2*X3|X2*X3)"}}
%---
%[output:5b67a257]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"(1)*(X1|X1)"}}
%---
%[output:0d57966f]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"(-1)*(X1|X1)"}}
%---
%[output:0fca7793]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"(1)*(1|X1*X2) + (1)*(X1|X2) + (-1)*(X2|X1) + (1)*(X1*X2|1)"}}
%---
%[output:12d298ff]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"(0)*(1|1)"}}
%---
%[output:6c0077b2]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"(1)*(1|X1*X2*X3) + (1)*(X1|X2*X3) + (-1)*(X2|X1*X3) + (1)*(X3|X1*X2) + (1)*(X1*X2|X3) + (-1)*(X1*X3|X2) + (1)*(X2*X3|X1) + (1)*(X1*X2*X3|1)"}}
%---
%[output:0f8e0a60]
%   data: {"dataType":"textualVariable","outputData":{"name":"Z","value":"    coeff    base\n    _____    ____\n      0       -  "}}
%---
%[output:1ed56b43]
%   data: {"dataType":"text","outputData":{"text":"X =\n  3 個のエントリをもつ <a href=\"matlab:helpPopup('dictionary')\" style=\"font-weight:bold\">dictionary<\/a> (<strong>double<\/strong> ⟼ <strong>VectExtAlg<\/strong>):\n    1 ⟼ 1×1 VectExtAlg\n    2 ⟼ 1×1 VectExtAlg\n    3 ⟼ 1×1 VectExtAlg\n","truncated":false}}
%---
%[output:53bfa2a8]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"    coeff    base\n    _____    ____\n      0       -  "}}
%---
%[output:81504f67]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"    coeff    base\n    _____    ____\n      1      x1x2"}}
%---
%[output:01c915eb]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"    coeff     base \n    _____    ______\n      1      x1x2x3"}}
%---
%[output:5a3c6adc]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"    coeff     base \n    _____    ______\n     -1      x1x2x3"}}
%---
%[output:9222b755]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"    coeff      base  \n    _____    ________\n      1      x1x2 ⊗ 1\n     -1      x2 ⊗ x1 \n      1      x1 ⊗ x2 \n      1      1 ⊗ x1x2"}}
%---
%[output:8083797b]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"    coeff       base   \n    _____    __________\n      1      x1x2x3 ⊗ 1\n      1      x2x3 ⊗ x1 \n     -1      x1x3 ⊗ x2 \n      1      x3 ⊗ x1x2 \n      1      x1x2 ⊗ x3 \n     -1      x2 ⊗ x1x3 \n      1      x1 ⊗ x2x3 \n      1      1 ⊗ x1x2x3"}}
%---
