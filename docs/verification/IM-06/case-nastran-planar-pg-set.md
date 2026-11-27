# IM-06 一つの Nastran 対応表にある全 Shell を準備する

認証済みの Nastran 対応表に含まれる面から、対象 Shell を初出順に選ぶ。同じ SI Case と認証デッキを渡して各 Shell の PG、厚さ、積層、材料軸を照合する。

すべて成功したときに各 Shell の所有結果をまとめて返す。[FD-03 対応表の契約](../../specs/foundation/external-mesh-sidecar-contract.md) [単一 Shell の準備](case-nastran-planar-pg.md)

## 所有結果と失敗時の再利用

一括入口で Case の親識別と物性を読んだデッキの出典を確認してから面を選ぶ。Shell の件数、選定時の面走査と既出 Shell との比較、Shell ごとの材料軸件数を限る。準備結果は単一 Shell の所有結果をそのまま保持する。

実デッキは Shell 2 の二つの三角形と、異なる PG を持つ Shell 1 の四辺形を含む。対応表では Shell 2、1、2 の順に面を宣言する。結果は Shell 2、1 の順で、それぞれの PG、元面、Case 層と共通のデッキ・対応表ダイジェストを確認した。

後で処理する Shell 1 の Case 軸だけを変えると、一括入口は Shell 1 とその宣言位置、元要素、Case 層を報告した。失敗後も同じ入力と作業領域で Shell 2 を準備でき、一括入口の再試行では同じ Shell 1 の不一致を返した。

選定は3面の走査と2回の既出比較を使い、計5回の上限で通過した。Shell 件数と選定回数が一件不足する場合も検査した。別 Case と事前取消は、共有 ResourcePlan に割当が増える前に拒否した。準備中に成功した割当は、後段で不一致を返した場合も共有 ResourcePlan に累積する。

## 検証に使った構成

対象ソースは `7996cbdab8ca98b1673bf7b6426969e387b403dd`。Windows 11 の MSVC 19.51.36260 Debug と Release、Ubuntu 24.04 の GCC 14.2.0 Debug（ASan・UBSan）で、`io.external_pg_regional_thickness`、`io.external_planar_pg`、`io.nastran_laminate_cards`、`io.reference_input` が各4件合格した。各構成で、実ファイルの対応表試験は1634アサーションだった。変更した翻訳単位の clang-tidy と clang-format 21.1.2 も通過した。
