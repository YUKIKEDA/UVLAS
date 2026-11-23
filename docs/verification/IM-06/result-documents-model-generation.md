# 失敗結果と確定状態のモデル世代

[結果出力の検証](result-documents.md)に加え、同じ確定IDを持つ異なるモデル世代の状態を照合した。失敗に記録した世代と一致する状態は受け付け、別世代の状態は `WriteResults` の `StaleGeneration` で拒否した。

Windows MSVC 19.51 Releaseで `io.result_documents` と `io.observation_csv` が通過した。結果JSONは5ケース122アサーション、観測CSVは5ケース120アサーション。観測CSVを含むPR差分のC++ファイルはclang-format 21.1.2の `--dry-run --Werror` を通過した。
