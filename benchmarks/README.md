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

## 密基準解法の容量baseline

`DensePreparation/{16,64}` は元係数・尺度のsnapshot、workspace、LU、条件推定を含む準備を測る。`DenseSolveBinary64/{16,64}` と `DenseSolvePromotedTwoWord/{16,64}` は準備済み分解を再利用し、元残差・後退誤差の認定と所有する出力の確保を含む。前者の対角成分は4、後者は3、右辺は全て1、尺度は1、各列は一つ。後者の許容差は1e-30、前者は1e-12であり、異なる精度要求の基準値として別々に読む。

いずれも合成小系の費用を追うbaselineで、実翼・実4096系・実務探索の性能受入値ではない。スレッド数は1、反復改善の上限は5、物理消費者checkは空。出力の確保はタイマーに含み、入力作成と最初の分解準備はsolveのタイマーに含めない。準備ケースは独立に測る。

`UpdatedPreparation/{16,64}` と `UpdatedSolve/{16,64}` は複素rank1更新の容量baselineである。$$ A_0=4I $$、$$ U $$ の全成分は $$ i/32 $$、$$ V $$ の全成分は $$ (1+i)/16 $$、右辺の全成分は $$ 1+i/8 $$、尺度は1、許容差は1e-12とする。prepareは元係数・U/V・尺度・workspace・実LU・複素小Schurを含み、solveは準備済み因子を再利用して元残差の認定と所有する出力の確保を含む。rank、実際の補正回数、出力確保の有無を各ケースのcounterへ記録する。

## 固定前処理GMRESの容量baseline

`gmres/Preparation/{96,1024}` と `gmres/Solve/{96,1024}` は、三重対角の合成系で準備と求解の費用を分けて測る。演算子は $$ A=\operatorname{tridiag}(-1/4,1,-1/4) $$、前処理は恒等作用、右辺と尺度は全て1とする。再始動幅は64、演算子作用の上限は512、右辺は1列、絶対残差と後退誤差の許容差は1e-11である。作用の丸め誤差の上界も元残差の認定へ含める。

準備は尺度のsnapshotとKrylov配列の確保を含む。solveは準備済み配列を再利用し、元残差・後退誤差の認定と所有する結果の確保を含む。入力作成と最初の準備はsolveのタイマー外である。物理観測の要求は空とし、実翼や実務規模の連成系の受入値には用いない。演算子作用数を各ケースのcounterへ残す。
