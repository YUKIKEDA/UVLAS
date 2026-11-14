# 0001: 依存はvcpkgのmanifestで固定する

- Status: Accepted
- Date: 2026-11-11
- Issue: #1

## Context

採用仕様N4-01は、Eigen 5.0.1、Spectra 1.2.0、Boost 1.92.0を「公式配布のURLとSHA256でロックし、未固定版の自動取得を呼ばず、ロック済みヘッダを指定する」と定めている。設計作業4では `build/numerical-deps/` に配布物と独自のロックファイルを置いて検証した。

製品の構築にこの独自形式を引き継ぐと、取得スクリプトとロック形式を自前で保守することになる。利用者は標準的な手法を求めた。

## Decision

依存はvcpkgのmanifest mode（`vcpkg.json`）で管理する。

- `builtin-baseline` でレジストリのコミットを固定し、`overrides` で各パッケージの版を固定する。
- Spectra 1.2.0はレジストリに無い（1.0.1まで）ため、`cmake/vcpkg-ports/spectra` のoverlay portで固定する。取得物のSHA-256は、設計作業4で使った配布物と同じ `fdfccacb…a73b`。
- 試験専用のCatch2と、Python結合のpybind11は、manifestのfeature（`tests`、`python`）に分ける。
- PresetがvcpkgのtoolchainをCMakeへ渡す。CMake側で依存を取得しない。

N4-01の「ロック済みヘッダを指定する」は、「vcpkgのbaselineとoverridesで固定した版を使う」と読み替える。版そのもの、ライセンスの扱い、最新版へ黙って追従しないという要求は変えない。

## Consequences

- 版の更新は `vcpkg.json` の変更として差分に現れる。更新には、数値結果と性能の再確認を伴う。
- 取得物の同一性は、vcpkgの各portが持つSHA-512で検査される。設計作業4の配布物との一致は、Spectraはハッシュで、他は版で確かめている。Eigen、Boostなどのバイト一致は確かめていない。
- 開発機とCIにvcpkgが要る。`VCPKG_ROOT` を設定する。
- 設計作業4の `build/numerical-deps/` とそのロックファイルは、当時の証拠として残す。製品の構築には使わない。
