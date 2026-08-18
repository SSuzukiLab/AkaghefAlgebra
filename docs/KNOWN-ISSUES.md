# 既知の破損と未整備

最終確認: 2026-08-18（MATLAB R2025b で実行確認）。

このファイルは**揮発性**。思想を書く他の docs とは性格が違うので分離してある。
修正したらここから消すこと。

## クラスがロード・実行できないもの

| 対象 | 症状 | 直し方 |
|---|---|---|
| `StrQUEAlg` および派生 (`StrAnAlg`, `StrCnAlg`) | `StrQUEAlg.m:20` の `rank`（Lie 環の rank）が `StrAlg.m:24` の Dependent `rank`（テンソル階数）と衝突し、クラス定義がロードできない | Lie 側をリネーム |
| `PolAlg` および派生 (`QWeylAlg`, `WeylAlg`, `PolCnmodalg` 等) | `PolAlg.m:468` が呼ぶ `sortrowCustom` が存在しない | `PolAlg.m:456-467` にコメントアウトされた旧実装がある |
| `StrWeylAlg` | `algbase` プロパティが未定義 | `base` に置換 |
| `StrAnAlg_` | 同上（`:20, :31`）、`set_cp` の引数過多（`:132`）、`dimV` 代入（`:133`） | 3 箇所修正で動作確認済み |

## 数学的なバグ

いずれも `Examples/StrAlgebra/StrQUEAlg.m` の `setRelation` 内。**未実行のため
発覚していない**（`checkRepresentation` の呼び出しはリポジトリ全体で 0 件）。

- `:163,165,167,169` — `dot(a(i),a(j))`。`a = CD.SimpleRoots` は
  $n\times\dim$ 行列なので `a(i)` は線形インデックスのスカラーになる。
  正しくは `dot(a(i,:),a(j,:))`。
- `:189` — `qNumS.binom(q^d(i),kMax,i)`。第 3 引数はループ変数 `k` であるべき
  ところが Dynkin ノード番号 `i` になっている。量子 Serre 関係式の係数が壊れる。

## 未整備

- **テストが実質ない。** `tests/HelloWorldTest.m` はスタブ、`tests/TestSuite.m` は
  `tests` を代入前に参照していて動かない、`tests/TestCalcSweedler.m` は unittest
  ではなく Live Script 形式のスクリプト。
  実質的な回帰テストは `Execution/C2410/Uqsl2Test.m`（69 行、assert ベース）と
  `Usl2Test.m` のみ。テストを整備するならこの 2 本が出発点。
- **設定が 2 系統に分裂している。** `AlgebraConfig.m`（`dynamicprops` シングルトン）と
  `CR.m`（別シングルトン）が役割を重複させている。`CalcRuleInit.m` は
  セッション開始時に手動実行が必要で、自動フックはない。
- **正本と実験用の実装が区別されていない。** 積公式がまだ見つかっていない段階で
  何通りか試した実験実装が、正本と並んだまま残っている（`QWeylAlg` の
  `prodevalQW1`/`prodevalQW2`、`PolCnmodalg`/`PolCnmodalg2` など）。
  実験を残すこと自体は妥当だが、どれが正本かを名前か配置で示す必要がある。
  詳細は [series-pol.md](series-pol.md) を参照。
- **antipode が全クラスで未実装。** `Str` 系の `S` は `% NG` コメント付き。
  `StrUqsl2.m:106-112` の `c=(-1)^length(p)` は $U_q$ では誤り。
- **`.mlx` のクラス名が追随していない。** `Execution/` 配下の C 型 Live Script が
  `symC2alg`, `symCnmodalg`, `strCnmodalg2`, `symp` を参照しているが、実体は
  `PolC2modalg`, `PolCnmodalg`, `PolCnmodalg2`, `StrCnmodalg`。
  `sym*` → `Pol*` のリネームに追随漏れ。Windows では大文字小文字非区別のため
  一部が偶然通る。
- `Core/qAnalog/qNumS.m:52-67` の `factor` は中身が全部コメントアウト。
- ルート直下の `C250814Uqsl2BorelSmallStr2Vec.mat` は 0 バイト。

## 動作確認が取れているもの

参考までに、${TODAY} 時点で動くことを実測した範囲。

- `Usl2`, `StrUqsl2`, `StrWeylXQ`, `StrCnmodalg`, `CartanData('A3')`, `qNumS`
- `StrUqsl2` の $q$ 微分表現に対する `checkRepresentation` — 定義関係式 7 本すべて成立
- `StrAnAlg_` を上表の 3 箇所修正した上で $A_2$（$\mathfrak{sl}_3$）に適用
  — 量子 Serre 関係式を含む 43 項目すべて成立、所要 32.5 秒

なお係数が 2 倍された偽の関係式では判定が落ちることも確認済みなので、
`checkRepresentation` の判定自体は機能している。ただし
[DESIGN.md](DESIGN.md) に書いたとおり、ゼロ判定は `simplify` 依存で決定手続きでは
ないことに注意。
