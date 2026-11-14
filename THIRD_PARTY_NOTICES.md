# 第三者ライセンス通知

UVLAS自身のソースのライセンスは、ルートの `LICENSE` にあるApache-2.0である。第三者のソース・ヘッダには、それぞれのライセンスが適用される。

| 依存 | 固定版 | ライセンス | 用途 |
|---|---|---|---|
| Eigen | 5.0.1 | MPL-2.0、および同梱の互換ライセンス | 数値計算。`EIGEN_MPL2_ONLY` を定義する |
| Spectra | 1.2.0 | MPL-2.0 | 一般固有値問題 |
| Boost Multiprecision／Hash2 | 1.92.0 | BSL-1.0 | 精密参照、有理数、SHA-256 |
| nlohmann/json | 3.12.0 | MIT | JSON |
| Catch2 | 3.16.1 | BSL-1.0 | 試験専用 |
| pybind11 | 3.1.0 | BSD-3-Clause | Python結合のfeature。現在のPresetでは導入しない |

上の表はライセンス全文や著作権表示の代わりではない。CMakeの構成時に、vcpkgが導入した各パッケージの `share/<package>/copyright` を、構成先の `third-party-notices/<package>.txt` へバイトを変えずに複製する。導入済みの推移的依存とビルド用パッケージの通知も含むため、この一覧は製品のリンク依存の一覧ではない。Eigenの通知には、同梱されたBSD等の全文も含まれる。

数値依存、JSON、試験を有効にしたときのCatch2、Python featureを有効にしたときのpybind11の通知が欠けていれば、構成を失敗させる。`toolchain.third_party_notices` が原本と複製のSHA-256一致を検査する。実際のCLI・wheel等への同梱と配布物検査は、ロードマップのIM06-0で実装する。

依存の取得・版の正本は `vcpkg.json` とADR 0001である。Spectraの取得物は設計作業4と同じハッシュで固定している。Eigen・Boost等は同じ版であるが、当時の配布物とのバイト一致は未確認であり、この機能試験によって同一と認定するものではない。
