# 数値基盤の凍結参照を製品経路へ取り込む

この文書はIssue #15の先行整備を定める。既存のPG／全自由後流24状態の参照を、元のバイト列と来歴を保持してfixtureへ取り込む。ここで確認するのは、固定された行列に対する製品の代数処理である。PG・UVLMの製品での組立、固有値・Arnoldiの製品経路やIM-00全体の完了を主張しない。

## 原本と正規の保存先

正本の索引は `rd/evidence-index.json`、引用元は `rd/foundation/numerics/records/numerical-foundation-adoption-audit.json` である。引用名 `build/analysis/numeric-oracles/` の実体は、`tools.rd.data_store.cached(name)` が返す `out/rd-evidence/` 内のcacheにある。作業ツリーのbuildに無いだけで欠落と判定せず、索引のbundleと正規cacheを確認する。

取り込む原本は `pg-free-step.txt`、`pg-free-oracle.npz`、`oracle.json` の3個である。索引・元auditの引用・cache実体のSHA-256を3者で照合してから読む。元のexport scriptとFreeProblemの生成ソースも、audit引用・現在file・fixture mapを3者で照合する。原本とR&Dの記録を書き換えない。

原本はfixtureのsource配下へバイト単位で複製する。`.gitattributes` はこの狭い範囲を `-text` とし、引用済みのCRLFやnpzを正規化しない。派生する `fixture.json` はUTF-8/LFとし、元recordのid/sha256/hint、元artifactのid/hash/論理cache名とbundle、map型のsources/artifactsを持つ。利用者名や実機の絶対pathは記録しない。

## 固定する代数問題と独立参照

元の行列Pの24×24係数は、textとnpzのbinary64 bitsが一致することを確認する。npzはこの原本のNPY-v1・little-endian binary64/complex128・row-major形式だけを読む。pickleや任意dtypeは扱わない。

製品求解へ渡す行列は、各成分のbinary64減算で定義した $$ A=I-P $$ である。rhsは $$ b_i=(i+1)/24 $$、状態の並びは抽象的な状態番号 0〜23、SI指数はdimensionlessの代数replayである。元の物理stateごとの未知の単位や尺度を推測して与えない。

原binary64係数をexactにDecimalへ取り込み、UVLAS/Eigen/Spectraを使わないscalarのpartial-pivot消去で参照解を120桁に固定する。160桁の別計算との差を $$ 10^{-100} $$ 以下と確認する。検証時は参照解を作り直さず、保存した参照の原残差を160桁で確認する。

製品試験はrow-major fixtureを列優先の所有DenseRealへ一度だけ変換し、DensePreparedとsolve_linearを実行する。相対解差と原残差は固定した $$ 10^{-12} $$ 級の代数予算で比較する。この予算は、凍結されたA/bに対するものであり、Pの物理的な微分・FD・格子誤差の精度ではない。npzの完全固有値は原本のまま保持し、和とtraceの整合を確認するが、製品固有値solverの成功へ数えない。

## 取込みと検証の入口

repo rootを作業directoryにして、`python -m tools.qa.foundation_reference import --output <new-directory>` を使う。cacheが別のcheckoutにある場合だけ `--cache-root` でcache rootを明示する。入力cacheは読取り専用で、既存output directoryは処理前に拒否する。新directoryをexclusiveに確保し、全原本の書込み後にfixture.jsonをcompletion markerとして作る。失敗時に既存の入力・出力を上書きしない。

`python -m tools.qa.foundation_reference verify <directory>` は、hash、id、shape、order、finite、JSON型、複素pair、固定予算、独立参照の残差を検査する。boolを数値0/1へ変換せず、全pairの実部・虚部を同じ規則で扱う。ソースとfixtureのhashだけを同時に付け替えても、元auditと一致しなければ拒否する。数値実験や物理モデルの再生成は行わない。

CTestも同じmodule入口を使い、WORKING_DIRECTORYをrepo rootへ固定する。source treeのpathが暗黙にsys.pathへ入る環境に依存しない。工具は標準libraryだけで動作し、製品や試験からrdのmoduleをimportしない。

## IM00-11の未達を残す

| 必須項目 | 現在の製品経路と、この先行PRでの状態 |
|---|---|
| 密実／複素更新・不定KKT・大規模疎 | IM00-6/7に製品実装と試験がある。このPRは新しい実PG由来の代数fixtureを密求解へ接続する。容量の既存試験は性能認定ではない |
| 一般実／複素の完全固有値、負の対称一般化固有値 | 製品入口・所有結果・原固有残差の認定が未実装。R&DとtoolchainのEigen実行は代わりにならない |
| 一般複素Arnoldiと実PG全履歴Pの固有照合 | Spectraのtoolchain試験はN4-01だけである。製品の作用素・所有結果・失敗／資源／世代契約と照合が必要で、凍結Pを輸入しただけでは達成しない |
| GMRES | Issue #12で実装中。完了をこのPRから推定しない |
| two-product／有理rank／区間平方根／SHA | EFTと状態fingerprintの製品処理は存在する。有理rankと256bit平方根囲いの既存試験は依存／R&Dの段階であり、製品consumerへの接続をgateで棚卸しする。このPRは新しい算術APIを追加しない |
| 型・世代・予算・投影・失敗後不変 | 既存製品試験を根拠ごとにgateへ集約する必要がある。新fixtureでは予算拒否後の入力・既存解の不変も確認する |
| 両OS clean build、4数値群の製品gate | 未達。fixtureのimportとWindows局所試験を、IM00-11完了へ置換しない |

固有値基盤の不足はIssue #15の中で明示する。N4-07と開始審査2節は小系完全スペクトルの照合と一般Arnoldiを要求している一方、現ロードマップのIM00-6/7/8にはその製品入口がない。後続IM-05の物理mode探索全体を前倒ししたり、未実装をtoolchain試験で埋めたりせず、製品gateを閉じる前に必要な基盤範囲を定める。

複素pairの型検査は全24組の実部・虚部に負対照を置く。最初のpairだけでは、その値が非零なため、boolの誤変換を型違反として検出したのか、単なる参照との差として検出したのかを区別できなかった。全組へ広げた試験は、boolを数値へ通す変異を8箇所で検出した。代表成分の検査を、全系列の型保証へ読み替えない。
