# 静的 Case の観測参照 JSON

R-05 の静的平衡で、点変位と支持反力を独立な参照値と比較する。Case に宣言した SHA-256 で、数値、単位、物理条件、比較条件を含む一つの UTF-8 JSON ファイルを照合する。この契約は [FD-09](case-schema-fields.md#fd-09-観測誤差配分) と [外部参照の取得契約](reference-snapshot-start-contract.md) に追加して適用する。

## FD-09 の参照形式

`uvlas.observation-reference/1` の初回適用は、[平面 PG の静的平衡](../structure/pg-flat-static-equilibrium-contract.md) における世界基底のシェル点変位と支持レンチとする。静的評価点では `instantaneous` と `right` を指定し、Case の比較窓の両端と評価時刻を一致させる。FD-09 の他の観測、時刻窓、基底にはこの版を適用せず、対応する参照形式を規定するまでファイルによる比較を `Unsupported` とする。

最初に重複キーのない JSON として構文を検査し、ルートの `schema` で版を判別する。未対応版はその版のキー構成を検査する前に `Unsupported` とする。この版の各オブジェクトは以下のキーだけを持ち、省略を許さない。既存の厳格な JSON 入力と同じ構文、重複キー検査、十進字句の桁数と指数の上限を適用する。`Number` は有限の JSON 十進数を厳密値として扱う。ID は [Case の共通型](case-schema-contract.md) の `Id<T>` と同じ十進文字列で記す。0 を含む uint64 の全域を扱う。

| オブジェクト | 必須キーと値 |
|---|---|
| ルート | `schema` は `uvlas.observation-reference/1`、`condition_sha256` は次節の小文字16進64文字、`observations` は1件以上の参照レコード |
| 参照レコード | `id` は Case の Observation ID、`observation` は下表の記述子、`comparison` は下記の評価条件、`components` は下表の名前から選んだ1件以上の成分 |
| 評価条件 | `statistic` は `instantaneous`、`time_s` は秒単位の Number、`time_side` は `right` |
| 成分 | `name` と `unit` は下表の文字列、`value` と `uncertainty` はその単位での Number。`uncertainty` は非負の絶対上限 |

| 観測 | `observation` のキーと値 | `components` の名前と SI 単位 |
|---|---|---|
| シェル点変位 | `kind` は `point`、`quantity` は `displacement`、`basis` は `{"kind":"world"}`、`point` は FD-02 の `MaterialPoint` と同じ形で `target.kind` を `shell` に限る | x、y、z。全て `m` |
| 支持レンチ | `kind` は `reaction`、`basis` は `{"kind":"world"}`、`support` は Support ID、`origin_m` は世界座標の3個の Number | Fx、Fy、Fz は `N`。Mx、My、Mz は `N*m` |

一つのファイルを複数の観測から参照できる。ファイル内に同じ ID は一度だけとし、このファイルを参照する Case の ID 集合に厳密に一致させる。レコードと成分の配列順は任意とする。Case の `ComparisonReference.columns` とファイル内の成分名はそれぞれ重複を拒否し、両者の集合を一致させる。列と不確かさは Case の列順で対応付け、比較結果は FD-09 の物理成分の順で返す。選択外の成分は未検証とし、R-05 の受入には表の全成分を要求する。点の対象・材料座標・厚さ側、または支持・モーメント原点を、各観測の宣言と照合する。原点と時刻を照合するときは、Case の長さと時間を SI へ厳密に換算する。点の材料座標は無次元のまま比較する。`time_side` は時刻の側を表し、厚さ側は `point.side` に保持する。

各成分の `ComparisonReference.uncertainty` を Case の根単位から SI へ厳密に換算し、ファイルの `uncertainty` に一致することを確かめる。読取りと照合は数値 kernel 用の丸めを行う前の厳密な値で実行する。

## Case の条件を照合する

`condition_sha256` は `schema`、`units`、`model`、`flow`、`analysis`、`numerics` と、存在する場合の `initial`、`restart` だけを持つ Case の投影を照合する。物理条件に関わる外部入力には元ファイルの宣言 digest が含まれ、参照取得時に実際のバイト列との一致を別に検査する。比較参照のパスと digest は投影に入らないため、参照 JSON を書いてから元 Case にその SHA-256 を記入できる。投影には宣言内容を使うため、モデルの幾何や材料が等価でも、内製、Gmsh、Nastran の各入力をそれぞれの digest に固定する。

ハッシュの入力は ASCII の `uvlas.observation-condition/1\n`（末尾は LF）に次の符号化を連結したバイト列とする。JSON 文字列は復号後の UTF-8 を使い、Unicode の正規化は行わない。長さと個数は先行零のない十進と ASCII `:` で表す。

- オブジェクトは `o`、キー数、キーと値の繰り返しの順に符号化する。キーは UTF-8 の符号なしバイト順で昇順に並べ、通常の文字列として符号化する。
- 配列は `a`、要素数、元の順の要素を連結する。文字列は `s`、長さ、UTF-8 バイト列を連結する。
- 数字は `d`、正規化した字句の長さ、字句を連結する。零の字句は `0`。それ以外は値を $c10^e$ とする符号付き整数係数 $c$ と指数 $e$ を 係数、小文字 `e`、指数の順で空白を挟まずに格納する。例えば `1.2500` は `125e-2` になる。$c$ の末尾に零を置かず、両整数は正符号と先行零を省く。
- 真、偽、null はそれぞれ ASCII 1バイトの `t`、`f`、`n` とする。

例えば符号化対象が `{"x":[1.25,"a"]}` のとき、版の接頭辞に続くバイト列は `o1:s1:xa2:d6:125e-2s1:a` となる。投影のハッシュは形式を固定する照合値であり、チェックポイントの互換性指紋とは別に管理する。

## 取得から確定観測の比較までをつなぐ

CLI のファイル入力と Python のファイル・辞書入力から、既存の Case 参照取得を呼ぶ。役割別の読取りは `observation-reference` として取得済みの元バイト列だけを消費する。閉じた構文を検査した後に依存辺が空であることを確定する。検査後のパス再取得を禁止する。宣言元と参照ファイル内の失敗位置の両方を診断に保持する。

参照検証に成功すると、親の `SchemaIdentity`、Source ID と SHA-256、condition SHA-256、観測 ID、照合した記述子・条件・各成分の厳密値を読取り専用で取得できる。結果は値を所有する。別の Case の結果と組み合わせた場合は `StaleGeneration` で拒否する。構文、digest、成分、条件の不一致は `InvalidInput` とする。未対応版や、この版の対象外の Case 観測は `Unsupported` とする。対応する Case 観測に対してファイルが別の側や基底を指定した場合は、入力の不一致として扱う。資源上限と取消は既存の失敗分類に従う。参照取得のファイル上限と既存 JSON のノード・深さ・文字列・数字の上限を適用する。ファイル内レコード数はそのファイルを参照する Case の観測数以下、成分数は点変位で3、支持レンチで6以下とする。同じファイルは一回だけ解析し、照合は全宣言に対して行う。元バイト列の解析、投影の走査・整列・ハッシュ、厳密数の換算、全参照の所有結果を checked 算術で見積もり、確保前に計上する。処理中と公開前に取消を確認する。失敗・取消の途中結果は公開しない。

数値比較では同じ Case と同じ準備済みモデルおよびレイアウトに属する確定状態から観測を評価する。観測結果の親 `SchemaIdentity` と準備済みモデル・レイアウトの同一性、確定状態の ID と世代を照合する。同じ数値の ID が別の所有者で使われた場合は拒否する。観測 ID、時刻、時刻側、基底・物理点または支持と原点を照合する。試行状態は受け付けない。成分 $i$ の製品値が全評価誤差を含む有限で順序付きの包絡 $[l_i,h_i]$、参照の厳密値が $r_i$、不確かさが $u_i$ のとき、誤差の上限は次式とする。

$$
E_i=\max\bigl(|l_i-r_i|, |h_i-r_i|\bigr)+u_i .
$$

各成分について、Case の SI 換算済み絶対許容値 $a_i$、相対許容値 $b_i$、正尺度 $s_i$ から作る $a_i+b_i s_i$ と、[A5-01](../../requirements/practical-acceptance-plan.md#a5-01-判定単位) が定める静変位・力・モーメントの固定尺度と許容値を両方適用する。Case の尺度を変更しても固定の上限は変わらない。誤差の上限を外向きに、許容値を内向きに囲うか、厳密な比較を使う。選択した全成分で上限が許容値以下と証明できた場合だけ、その成分集合を合格とする。観測が得られない場合、包絡が成立しない場合、誤差が予算を満たさない場合をそれぞれ記録する。参照検証の結果と、求解後の観測比較の合否を別に報告する。FD-09 の適用誤差群と全数値予算の認定も別に行う。参照の不確かさは参照比較に一度だけ計上する。

## 矩形 Case の参照ファイル例

手入力の矩形 Case で変位と支持レンチを同じファイルから参照する例を示す。Case の二つの観測が同じ Source ID、パス、元バイト列の SHA-256 を宣言する。Case の `columns` はそれぞれ `["x","y","z"]` と `["Fx","Fy","Fz","Mx","My","Mz"]` とし、対応する不確かさは全て零とする。下記の整形を変えたときは、元バイト列の SHA-256 を再計算する。

```json
{
  "schema": "uvlas.observation-reference/1",
  "condition_sha256": "378fea09ac121f7e2240224bfacb63c983ab73b695fe39199ab710685c900b0e",
  "observations": [
    {
      "id": "1",
      "observation": {
        "kind": "point",
        "quantity": "displacement",
        "basis": {"kind": "world"},
        "point": {
          "target": {"kind": "shell", "id": "1"},
          "coordinates": [1, 0.5, 0],
          "side": "interior"
        }
      },
      "comparison": {"statistic": "instantaneous", "time_s": 0, "time_side": "right"},
      "components": [
        {"name": "x", "unit": "m", "value": 0, "uncertainty": 0},
        {"name": "y", "unit": "m", "value": 0.001953125, "uncertainty": 0},
        {"name": "z", "unit": "m", "value": 0, "uncertainty": 0}
      ]
    },
    {
      "id": "2",
      "observation": {
        "kind": "reaction",
        "basis": {"kind": "world"},
        "support": "1",
        "origin_m": [2, 1, 3]
      },
      "comparison": {"statistic": "instantaneous", "time_s": 0, "time_side": "right"},
      "components": [
        {"name": "Fx", "unit": "N", "value": 0, "uncertainty": 0},
        {"name": "Fy", "unit": "N", "value": -993641.3727700710296630859375, "uncertainty": 0},
        {"name": "Fz", "unit": "N", "value": 0, "uncertainty": 0},
        {"name": "Mx", "unit": "N*m", "value": 993641.3727700710296630859375, "uncertainty": 0},
        {"name": "My", "unit": "N*m", "value": 0, "uncertainty": 0},
        {"name": "Mz", "unit": "N*m", "value": 993641.3727700710296630859375, "uncertainty": 0}
      ]
    }
  ]
}
```

## 実装の検証

- R-05 の手入力、Gmsh、Nastran の三経路から参照元と Case を独立に生成する。矩形の解析解から、先端変位 $1/512\ \mathrm m$ と、世界座標 $(2,1,3)\ \mathrm m$ を基準点とする支持レンチを得る。$F_y=-4167633984375/4194304\ \mathrm N$、$M_x=M_z=4167633984375/4194304\ \mathrm{N\,m}$。他の成分は零。解析値の十進字句・成分順・単位と独立 SHA-256 を照合する。条件の投影・符号化も製品から独立に計算する。
- CLI と Python のファイル入口で三経路の参照を照合する。Python の辞書入口にも明示した参照根を渡し、同じ条件ハッシュになる例には二進数で表現できる辞書値を使う。厳密値が異なる辞書入力には、実際の辞書値に基づく条件ハッシュで別に照合する。キー順、無関係な空白、等価な十進字句を変えたときの条件 digest の一致も試す。
- 参照値のファイルを変更したときに、元の digest では失敗することを確かめる。異なる Case 条件、点または支持、原点、基底、物理単位、時刻・時刻側、選択成分の欠落・重複、および列と不確かさの対応、負の不確かさ、未知の版をそれぞれ変更し、再計算したファイル digest を与えて契約違反の位置を照合する。
- 観測比較では参照と一致する確定値、許容値と等しい上限、許容値を超える上限、不確かさだけで許容値を超える場合、ID が同じ別の親や準備済みモデル、誤った確定世代を検査する。静的求解の受入では、元残差・拘束とともに、解析値との全成分の比較を実行する。資源上限、取消、失敗後の親の再利用も検査する。
