# Algebra

Public repository of MATLAB tool for computing mathematical objects.

生成元と関係式で定義された非可換代数、$q$ 微分作用素、有限次元 Hopf 代数などを
MATLAB (Symbolic Math Toolbox) 上で計算するためのライブラリ。

## はじめに読むもの

このライブラリには `Str` / `Pol` / `Vec` という三つの系列があり、どれも「代数の元」を
表現する。**まず [docs/DESIGN.md](docs/DESIGN.md) を読むこと。** 三系列がなぜ分かれて
いるか、どれを選ぶべきかが書いてある。ここを知らずにコードから入ると設計を読み違える。

| ドキュメント | 内容 |
|---|---|
| [docs/DESIGN.md](docs/DESIGN.md) | 全体の設計思想、三系列の選択基準、共通の約束事 |
| [docs/series-str.md](docs/series-str.md) | `Str` — 語と書き換え。生成元と関係式から始めるとき |
| [docs/series-pol.md](docs/series-pol.md) | `Pol` — 指数ベクトルと積公式。速度と記号指数が要るとき |
| [docs/series-vec.md](docs/series-vec.md) | `Vec` — 構造定数。有限次元 Hopf 代数を扱うとき |
| [docs/KNOWN-ISSUES.md](docs/KNOWN-ISSUES.md) | 現在動かない箇所（揮発性） |

使い方のチュートリアルと内部アルゴリズムの解説は既存の文書にある。

- `Execution/Docs/StrAlgInstructions.md` — $U(\mathfrak{sl}_2)$ を例にした `Str` の操作
- `Execution/Docs/W241223QiitaMATLABAdvent.md` — 書き換えアルゴリズムと MATLAB 実装技法

## 要件

- MATLAB R2023a 以降（`dictionary`, `combinations`, `arguments` ブロックを使用）
- **Symbolic Math Toolbox（必須）**
- 他のツールボックスは不要

セッション開始時に `CalcRuleInit` を実行して表示設定を初期化する。

## 構成

```
Core/
  StrAlgebra/   Str 系列の中核（StrAlg, StrEndV, CartanData, HopfAlg）
  PolAlgs/      Pol 系列の中核（PolAlg）
  VectAlgebra/  Vec 系列（双対代数、Drinfeld/Heisenberg double）
  common/       共通基盤（VectAlg, 抽象インターフェース, 型パラメータ）
  Base/         基底と数値型
  qAnalog/      q 数、q 二項係数、q 指数
  sequence/     二項係数、Stirling 数ほか
  tensor/       疎テンソルと Sweedler 記法（別系統）
Examples/       具体的な代数の実装
Execution/      計算ノート（Live Script）と文書
tests/          未整備。KNOWN-ISSUES.md 参照
```
