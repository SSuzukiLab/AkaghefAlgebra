
# StrAlg説明書

 $U(sl_2 )=\mathbb{R}[E,F,H]/Relation$ の具体例に沿って使い方の説明を行う

 $$ Relation=([H,E]-2E,[H,F]+2F,[E,F]-H) $$ 

便宜上の生成元の番号: $(1,2,3)\leftrightarrow (E,F,H)$ 


この番号について多項式が昇順になるように整列して計算,出力する。

# Define

生成元を定義する

```matlab
[O,E,F,H]=Usl2.getGenerator;
```

```matlabTextOutput
警告: クラス 'sym' はクラス フォルダー内で定義されているため、MATLAB パス上で前の位置にある同じ名前の関数よりも優先されます。将来のリリースでは、クラス 'sym' は優先されなくなります。 

競合する項目の場所については、こちらをクリックしてください。
この警告を回避するためのガイドラインについては、こちらをクリックしてください。
```
# Construct

係数と生成元の番号を入力する

```matlab
A=O.make([3 2],{[] [1]})
```

```matlabTextOutput
A = 
    coeff    base
    _____    ____
      3       1  
      2       E  
```

```matlab
B=O.make(3:4,{[] [2 3]})
```

```matlabTextOutput
B = 
    coeff    base
    _____    ____
      3      1   
      4      F H 
```
# add
```matlab
A+B
```

```matlabTextOutput
ans = 
    coeff    base
    _____    ____
      6      1   
      2      E   
      4      F H 
```

```matlab
C=A-B
```

```matlabTextOutput
C = 
    coeff    base
    _____    ____
      2      E   
     -4      F H 
```

```matlab
A-A
```

```matlabTextOutput
ans = 
    coeff    base
    _____    ____
      0       1  
```

```matlab
O+A
```

```matlabTextOutput
ans = 
    coeff    base
    _____    ____
      3       1  
      2       E  
```

```matlab
A+3
```

```matlabTextOutput
ans = 
    coeff    base
    _____    ____
      6       1  
      2       E  
```

```matlab
4+A
```

```matlabTextOutput
ans = 
    coeff    base
    _____    ____
      7       1  
      2       E  
```
# multiply

lb(x,y)=\[x,y\]はLie bracket

```matlab
F*E*2
```

```matlabTextOutput
ans = -2*H +2*E*F
```

```matlab
assert(H^2*F-4*F+4*F*H-F*H*H==O)
F^3*H^2*E^2
```

```matlabTextOutput
ans = -96*F*H +96*E*F*F +48*F*H*H -48*E*F*F*H +42*F*H*H*H +16*E*E*F*F*F -42*E*F*F*H*H +6*F*H*H*H*H +8*E*E*F*F*F*H -6*E*F*F*H*H*H +1*E*E*F*F*F*H*H
```

カシミール作用素 $EF+FE+\frac{1}{2}H^2$ は中心に属する

```matlab
lb(E^3*F*H*F*E^3,F*E+E*F+H^2*(1/2))
```

```matlabTextOutput
ans = +0*1
```
# TensorProduct

テンソル積は縦棒 | で記述する。


テンソル積 | より+,\-の演算子のほうが優先度が高いので注意する。

```matlab
-(H|F)+(F*E|F)
```

```matlabTextOutput
ans = 
    coeff     base  
    _____    _______
     -2      H ⊗ F  
      1      E F ⊗ F
```

```matlab
(F|F)*(E|E)
```

```matlabTextOutput
ans = 
    coeff      base   
    _____    _________
      1      H ⊗ H    
     -1      H ⊗ E F  
     -1      E F ⊗ H  
      1      E F ⊗ E F
```

```matlab
lb(H|H,E|F)
```

```matlabTextOutput
ans = 
    coeff     base  
    _____    _______
     -4      E ⊗ F  
      2      E ⊗ F H
     -2      E H ⊗ F
```

```matlab
(1|E)+(E|1)
```

```matlabTextOutput
ans = 
    coeff    base 
    _____    _____
      1      1 ⊗ E
      1      E ⊗ 1
```

```matlab
1|E+E|1 % 1⊗(E+E⊗1)と解釈されてエラーの可能性もある
```

```matlabTextOutput
ans = 
    coeff        base     
    _____    _____________
      2      1 ⊗ E ⊗ 1 ⊗ 1
```
# comutiply
 $$ \Delta (EF)-\Delta (E)\Delta (F)=0 $$ 
```matlab
Delta(E*F)-Delta(E)*Delta(F)
```

```matlabTextOutput
ans = 
    coeff    base 
    _____    _____
      0      1 ⊗ 1
```

```matlab
Delta(E^3)
```

```matlabTextOutput
ans = 
    coeff      base   
    _____    _________
      1      1 ⊗ E E E
      3      E ⊗ E E  
      3      E E ⊗ E  
      1      E E E ⊗ 1
```

```matlab
[O,E,F,K,Ki]=Uqsl2.getGenerator;
```

```matlabTextOutput
クラス 'Uqsl2' のメソッド、プロパティまたはフィールド 'make_' が認識されません。
エラー: Uqsl2.make (行 38)
            obj=Uqsl2().make_(cf,pw,Uqsl2.B);
エラー: Uqsl2.getGenerator (行 41)
            O=Uqsl2.make(0,{});
```

```matlab
F*E*F
```
