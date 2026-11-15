# kernelの小規模計測

最初の対象は製品のchecked_extent_productである。Google Benchmark 1.9.5を使い、結果はBenchmarkのJSONとして保存する。時間値をCIの合否や実務受入の時間制限には使わない。CIは変更の影響があるときだけLinuxの専用featureを構成・コンパイルし、kernelの一覧を確認する。

## Windows

通常の試験構成と異なるprefixへ導入する。既にある通常prefixのfeature集合を切り替えない。

```powershell
pwsh tools/dev/msvc.ps1 cmake -S . -B out/build/windows-bench -G Ninja -D CMAKE_BUILD_TYPE=Release -D CMAKE_TOOLCHAIN_FILE="$env:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" -D VCPKG_MANIFEST_FEATURES=benchmarks -D UVLAS_BUILD_BENCHMARKS=ON -D BUILD_TESTING=OFF -D VCPKG_INSTALLED_DIR="$PWD/out/vcpkg_bench_installed"
pwsh tools/dev/msvc.ps1 cmake --build out/build/windows-bench --target uvlas_kernel_bench
$env:UVLAS_BENCH_MACHINE = 'bench-host-a'
$env:UVLAS_BENCH_CPU = (Get-CimInstance Win32_Processor | Select-Object -First 1).Name
$env:UVLAS_BENCH_MEMORY_BYTES = [string](Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory
.\out\build\windows-bench\benchmarks\uvlas_kernel_bench.exe --benchmark_repetitions=3 --benchmark_out=out/kernel-bench.json --benchmark_out_format=json
```

各コマンドの終了コードを確認してから次へ進む。機械ラベルは個人名やhostnameを避け、CPUとメモリは実際の値を使う。

## Linux

```sh
cmake -S . -B out/build/linux-bench -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE="$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" -DVCPKG_MANIFEST_FEATURES=benchmarks -DUVLAS_BUILD_BENCHMARKS=ON -DBUILD_TESTING=OFF -DVCPKG_INSTALLED_DIR="$PWD/out/vcpkg_bench_installed"
cmake --build out/build/linux-bench --target uvlas_kernel_bench
export UVLAS_BENCH_MACHINE=bench-host-a
export UVLAS_BENCH_CPU="$(LC_ALL=C lscpu | sed -n 's/^Model name:[[:space:]]*//p' | head -1)"
export UVLAS_BENCH_MEMORY_BYTES="$(( $(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) ))"
./out/build/linux-bench/benchmarks/uvlas_kernel_bench --benchmark_repetitions=3 --benchmark_out=out/kernel-bench.json --benchmark_out_format=json
```

Linuxの採取手順はこの開発機では未実行。既存の環境にあるGCC 14・vcpkgを用い、道具を無断で追加しない。

## 読む順序

まずcontextのcommit、dirty、構成、命令セット、機械条件、設定とworker数を確認する。同じ条件の反復のばらつきを見てから比較する。プロファイラの手順は [kernel-profile Skill](../.agents/skills/kernel-profile/SKILL.md) を使う。測った時間を改善の証拠にする場合は、変更前後の結果の同一性も必要になる。

Benchmarkの生JSONにはexecutableの絶対パス等が入り得る。生データはout/に置き、公開する要約は別ファイルに必要項目だけを抽出して作る。既存の証拠を伏せ字のために上書きしない。
