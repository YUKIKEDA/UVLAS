---
name: kernel-profile
description: "UVLASの実装済みkernelをWindows WPR/WPAまたはLinux perfで採取し、支配する呼出しと計測条件を読む。対象はbenchmarksの実kernel。最適化の採否や未実装のソルバー性能の予測には使わない。"
---

# kernelのプロファイルを採る

最初の対象はbenchmarks/dimensions_bench.cppから呼ぶ製品のchecked_extent_product。ビルドと機械条件の設定は [benchmarks/README.md](../../../benchmarks/README.md) に従う。性能改善を採用する場合はsolver-performanceの比較と同値性確認を併用する。

## 採取の前に

1. cleanなcommit、最適化構成、命令セット、worker数、機械条件、入力を固定する。採取用の構成も固定する。RelWithDebInfoとReleaseは同じ最適化条件とは限らず、採取用構成の時間をReleaseの性能値へ読み替えない。
2. 対象を1回だけ短く起動し、kernelが走ることを確認する。道具の存在だけを採取成功と書かない。
3. 他のビルドやベンチを止める。WPRの既存sessionがあれば所有者を確認し、別の採取を勝手に停止しない。
4. 出力は新しいout/配下の名前にし、既存の採用証拠・計測記録を上書きしない。

## Windows：CPUの呼出しstackを読む

WPR/WPAの存在と `wpr -profiles` のCPU profileを確認する。必要な権限が無ければ不足として報告し、権限変更や道具の追加を無断で行わない。採取用にはREADMEと同じ構成引数で、ビルド先をout/build/windows-profile、CMAKE_BUILD_TYPEをRelWithDebInfoにしてビルドする。下の実行ファイルもその構成を使う。

```powershell
wpr -status
wpr -start CPU -filemode
.\out\build\windows-profile\benchmarks\uvlas_kernel_bench.exe --benchmark_filter=checked_extent_product --benchmark_min_time=10s
wpr -stop out/kernel-cpu.etl
wpa out/kernel-cpu.etl
```

各終了コードを確認する。開始できなかった採取に対してstopを呼ばない。開始した採取は対象が失敗しても、そのsessionを停止して記録を閉じる。

WPAでは対象processに絞り、CPU Usage (Sampled)のstackを表示する。checked_extent_product自身のサンプルと呼出し元を含むサンプルを分ける。準備・Benchmark集計・OS処理が支配していれば、kernelが支配すると結論しない。未解決symbolが多ければPDBの対応を確かめて採り直す。

## Linux：perfの呼出しstackを読む

既存のperfが使え、kernelの権限設定が採取を許すことを確認する。権限不足は報告し、perf_event_paranoid等を無断で変更しない。READMEと同じ構成引数で、ビルド先をout/build/linux-profile、CMAKE_BUILD_TYPEをRelWithDebInfoとしてビルドする。下のコマンドはそのsymbol付き実行ファイルをDWARF stackで採取する。

```sh
perf record --call-graph dwarf -o out/kernel-perf.data -- ./out/build/linux-profile/benchmarks/uvlas_kernel_bench --benchmark_filter=checked_extent_product --benchmark_min_time=10s
perf report --stdio -i out/kernel-perf.data
```

自己時間と子の呼出しを含む時間を分け、サンプル数と未解決stackを確認する。採取できなかったものを0サンプルの成功にしない。オプションの根拠は [perf record](https://man7.org/linux/man-pages/man1/perf-record.1.html) と [perf report](https://man7.org/linux/man-pages/man1/perf-report.1.html) のmanual。

## 返すもの

実行コマンドと条件、採取の成否、支配する関数と自己／包含の区別、未解決symbol、次に確かめる仮説を報告する。この小kernelから全ソルバーの時間を推測しない。

生のETL・perf.data・Benchmark JSONにはパスや他processの情報が含まれ得る。公開前に必要項目だけの別要約を作る。大きな証拠の保存先は所有者に確認し、採用済みの生データを後から書き換えない。

## 今回確認した範囲

2026-11-15、WindowsでWPR/WPAの存在と `wpr -profiles` のCPU profileを確認した。WPRの採取・WPA解析・Linux perfは未実行。手順の記載を実行済みの証拠には数えない。