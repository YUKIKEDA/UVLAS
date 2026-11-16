# IM-00から物理・入出力パッケージへの引継ぎ審査

## 基盤の再利用は可能、次の実装は子作業の契約を閉じてから始める

2026-11-16。審査対象は `1a3c98f` に採用された製品コード、試験、追補を含む検証記録である。Issue [15](https://github.com/YUKIKEDA/UVLAS/issues/15)・[16](https://github.com/YUKIKEDA/UVLAS/issues/16) の終了判定のうち、ここでは型・状態・所有・精度の引継ぎを扱う。数値系列の総合審査と最終CIは、同Issueの終了記録に対応づける。この審査では試験を再実行していない。

IM-00は、宣言された入力・尺度・世代・資源・精度政策の範囲で、後続の実装が利用できる共通基盤を提供している。ただし終了審査で、precision契約とdense契約がIM00-11へ残した実4096全系の製品認定が未了と判明した。既存row witnessは全系の投影拒否を示すが、全4096行の保持解の認定を代替しない。この未達を閉じるまでIM-00終了と後続パッケージの実装開始を判定しない。以下は既存基盤の再利用範囲を示すもので、契約の具体化は並列に進められる。後続の物理残差、観測値、幾何、再開データをこの基盤へ結び付ける責任は、それぞれのパッケージに残る。

この判断は「実装が採用仕様に適合する」という範囲である。内部のErrorCertificateや解析例の合格を、全用途の数値認定、実験妥当性、実務規模の性能、実務受入へ読み替えない。ロードマップの現在行と状態は [roadmap.md](../roadmap.md) だけが持つ。

## 製品試験が守る境界と、後続が引き受ける責任

根拠は [N4](../specs/foundation/numerical-foundation-adoption-specification.md)、[保持契約](../specs/foundation/numerical-storage-contract.md)、[状態契約](../specs/foundation/state-contract.md)、[Python結合契約](../specs/foundation/python-binding-contract.md) である。表の試験は代表的な所在で、各記録の凍結版と後続の訂正を併せて読む。

| 境界 | 製品で確認した内容と根拠 | 後続で必要なこと |
|---|---|---|
| 強い型・形状 | IDの暗黙混同を拒否し、int64寸法、積、byte数、外部int32境界を検査する。密・CSC・点・自由度情報は所有する。[checked_dimensions_test.cpp](../../tests/core/checked_dimensions_test.cpp)、[storage_test.cpp](../../tests/numerics/storage_test.cpp)、[im00-2](../verification/IM-00/im00-2.json)・[im00-4](../verification/IM-00/im00-4.json) | DofLayoutの添字選択は非線形拘束写像ではない。物理量の成分順、単位、基底、対応と写像を各kernelの契約で定める。 |
| 所有・世代 | Model定義、試行状態、評価対象、資格packetを所有し、別Model・別親・別世代、変更後に同じbyteへ戻した試行を拒否する。全cache項目を照合する。[state_test.cpp](../../tests/core/state_test.cpp)、[linear_validation_test.cpp](../../tests/numerics/linear_validation_test.cpp)、[im00-9](../verification/IM-00/im00-9.json) | 状態payloadはopaqueなbyte列である。q、速度、後流、履歴、物理cacheの内容と更新規則は未定義のまま流用しない。 |
| 資源・失敗 | 確保前のchecked計上、入れ子予約、失敗時の自身の予約解放と内側成功の計上保持、typed Failureと内側詳細の保持を試す。[resource_plan_test.cpp](../../tests/core/resource_plan_test.cpp)、[failure_test.cpp](../../tests/core/failure_test.cpp)、[im00-2-reentrancy](../verification/IM-00/im00-2-reentrancy.json)・[im00-3](../verification/IM-00/im00-3.json) | ResourcePlanは累積payload計上であり、外部ライブラリの全確保やRSSのhard capではない。実caseのgraph・fill・準備費用と外部backendの予算を別に検証する。 |
| 保持・投影・消費先 | TwoWord両語の原残差、実返却配列の再認定、必要consumer全件の文脈・予算・達成状態を照合する。全体荷重の成功で別の局所荷重を代用しない。[precision_test.cpp](../../tests/numerics/precision_test.cpp)、[linear_result_test.cpp](../../tests/numerics/linear_result_test.cpp)、[dense_contract_test.cpp](../../tests/numerics/dense_contract_test.cpp)、[im00-5](../verification/IM-00/im00-5.json)・[im00-77](../verification/IM-00/im00-77.json) | 任意の物理kernelが両語を消費することまでは保証しない。応力、仕事、後流、再開等の各consumerが実値から独自に認定する。 |
| 親と既存結果の不変 | 資格不足、別親、取消、資源不足、callback失敗、後続solveでも親・既存結果を変更しない。確定は一回の所有移転である。[state_test.cpp](../../tests/core/state_test.cpp)、[dense_contract_test.cpp](../../tests/numerics/dense_contract_test.cpp)、[im00-9](../verification/IM-00/im00-9.json) | 状態確定の共通機構に、物理的な残差・幾何・時間・精度の全必要判定を接続する。状態機構自体が物理適格性を判定するわけではない。 |
| workspaceと並列利用 | 同じWorkspaceの再入・同時可変共有を拒否し、例外でもleaseを解放する。返値は再利用scratchを借用しない。[state_test.cpp](../../tests/core/state_test.cpp)、各求解器の契約試験、[im00-9](../verification/IM-00/im00-9.json) | 独立workspaceの利用と、全物理処理の並列安全性・決定性は別に検証する。共有するprovider、cache、親の不変条件を子作業に記す。 |
| Python寿命・可変性 | 入力snapshot、readonly所有view、派生view・memoryview・GC、書込み再有効化の拒否、mutable操作の排他、成功包装前の準備を確認する。[test_ownership.py](../../tests/python/test_ownership.py)、[im00-10-ci-acceptance](../verification/IM-00/im00-10-ci-acceptance.json) | 対象は既存の保持・状態処理、通常GILの単一interpreterである。物理solve、case辞書、run/resume、wheel導入の実装・認定は含まない。 |

有理rankの微小差と256bit平方根囲いは [im00-77](../verification/IM-00/im00-77.json) で既存製品の保持・消費・拒否・包絡へ接続した。これは任意有理行列のrank APIや256bit製品平方根APIを追加したという意味ではなく、PGの実拘束rankを代替しない。

## IM-01・IM-03・IM-06へ渡せるもの

依存関係は [開始審査](implementation-start-review.md) と各Issueに従う。次の表は子作業を分解するための候補であり、新しいロードマップ行や実装開始の決定ではない。

| パッケージ | 利用できるIM-00基盤 | 最初の子作業で閉じる契約とgate |
|---|---|---|
| [IM-01 PG参照係数・疎拘束](https://github.com/YUKIKEDA/UVLAS/issues/17) | checked寸法、単位・尺度付きlayout、所有配列、精度保持、疎解法、世代・資源・失敗の共通機構 | 採用PGの参照係数と正確rankの責任範囲、位置・厚さ・層・面のtrace、局所持上げ、疎chartの入出力を定める。実接続graphによるrank・fill・準備容量を早期gateに置く。微小方向を落としてrankを合わせない。 |
| [IM-03 UVLM](https://github.com/YUKIKEDA/UVLAS/issues/19) | 密保持・基準解法、作用とGMRES、TwoWord、共通認定、workspace、工程計測 | 先行候補は有限渦kernelと共有辺の幾何である。採用式の向き・循環符号・コア・近接/自己作用、原作用とJVP/VJP、荷重の仕事、誤差と資源の契約を閉じる。構造との写像はIM-01/02に依存し、後流・全履歴の所有と費用は後続の明示した子作業で検証する。 |
| [IM-06 入出力・利用・配布](https://github.com/YUKIKEDA/UVLAS/issues/22) | 所有Model/状態、Failure、ID・hash、既存Python結合、ビルド情報と通知 | 先行候補は採用case形式の構文・参照検査と、CLI/wheelの小さい実ランタイム検査である。物理Modelへの変換、run/resume、結果・checkpointは対応する物理契約と結ぶ。開発機でのimport成功を配布互換性としない。 |

IM-06のschema名、必須ブロック、厳格な型検査、単位と相対参照、CLIの役割、Python対象、読込み対象は [入出力採用仕様](../specs/foundation/io-execution-adoption-specification.md) に決定済みである。初期子作業ではその全体と矛盾しない検査範囲を定める。物理内容を欠いたModelのbyte列だけで、caseの意味検査や再開互換を完了したとしない。

IM-01の参照・rank契約、IM-03の有限渦kernel契約、IM-06の入力/配布基盤契約は、共通型の再設計を待たず並列に具体化できる。構造空力写像や物理入力変換には前段の契約を渡す。未確定の依存を仮の物理モデルで埋めない。

## Grillと、実装者が閉じる契約を分ける

既決の方式・精度・初回実務範囲を再選定する必要はない。C++23、採用シェル、採用UVLM、全必要観測と精度政策は維持する。まず各Issueを、依存が閉じる子作業とその検証単位へ分解する。

利用者へのGrillは、採用文から一意に決まらず、利用手順・対応範囲・受入条件を左右する所有者の決定に限る。たとえばREL-1（Issue #26）の配布対応環境・公開先・署名が該当する。独立作業の優先順は、依存関係と要求への効果から実装担当が決め、並列に進める。初回必須範囲は削らない。質問が必要な場合は、相互依存しない事項をまとめ、各選択肢・推奨案・メリット・デメリットを提示する。

式から導ける成分順やchecked容量、内部の型名・配置・routine分割は実装担当が根拠を示して決める。各子作業には、少なくとも次を明文化する。

- 入力から要求観測までの式、離散残差、JVP/VJPまたはtrace/仕事の関係。原問題と近似問題を区別する。
- shape、成分順、単位、基底、時刻、向き、所有と借用、親・試行・cacheの世代、変更後の無効化。
- 保持精度、原残差と観測予算、必要consumer全件、明示投影の再認定、失敗分類と親不変。
- checked準備容量、外部確保の保証範囲、工程計測、独立oracle、負対照、意味ある変異、実graphまたは実kernelで止めるgate。

以上をIssueの完了条件と試験へ対応づけ、Grill後にロードマップの行を追加してから実装へ進む。IM-00の終了だけを根拠に、未分解のIM-01・IM-03・IM-06を一括実装しない。
