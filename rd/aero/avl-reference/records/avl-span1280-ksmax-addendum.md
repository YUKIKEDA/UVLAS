# N1280格子生成の局所容量：追加事前仕様

2026-10-31。最初の [N1280実行記録](avl-span1280-fine-runs.json) は、FW-V01・M3のLOAD段階で `MAKESURF: Array overflow. Increase KSMAX to 640` により停止した。Fortran STOPの終了コードは0だったが、console検査で拒否した。解・荷重出力はなく、正常結果に数えない。元の実行・準備・ソースとプロトコルは保持する。

前回のNVX／NVMAX／NSMAXだけの調査は不十分だった。amake.fの局所KSMAX=500はNSMAXから独立している。YPT・YCP・YZLEN・IPTLOCと派生KPMAXがこの値を使う。今回は別コピー `source-symmetry16000-ks1000-spacing-pi64` で**amake.fのKSMAX=500→1000だけ**を変更する。KCMAX=50と数値処理は維持し、他の局所容量もソース監査する。COMMON／moduleを変更しないのでamake.oのみ再コンパイルし、残り40オブジェクトとライブラリは親版とハッシュ一致を要求して再リンクする。

新たな8ケース回帰は [容量仕様](avl-span1280-capacity-experiment.md) の格子・設定・時間・メモリ・保存量・バイト一致条件を維持する。比較対象は親16000版の今回の8回帰。回帰前に8つのN1280デッキをLOADしてQUITする格子生成試験を行う。解く操作Xは含めない。各10秒・合計60秒・RSS256MiB・出力8MiB以下とし、consoleの未知数・帯数が入力と一致し、overflow／STOP由来の欠落を検出する。最後のQUITメニューへ戻ったことを確認する。小弦だけの読込みで最細弦を代替しない。

回帰・格子生成・ソース監査が成立した後のN1280は別ランナー／別JSONを使用する。[細分化事前仕様](avl-span1280-fine-experiment.md) の順序・8ケース・資源・不変性・保存観測・収束条件は変更しない。スパンのN160／320／640には今回のKSMAX修正版の回帰データを使用し、同じ版で揃える。初回失敗を上書きしない。
