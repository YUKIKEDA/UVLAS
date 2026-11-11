# Linux C++23検証環境

通常のGCC／Clangと標準ライブラリが必要機能を満たす環境では、専用のパッケージ管理ツールは不要。以下はUbuntu 22.04のGCC 11を変更せず、C++23を検証した環境の再現手順。

使用ツールはmicromamba 2.10.0。取得方法は [公式手順](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) を参照。シェル初期化やシステムへのインストールは必要ない。

`linux-gcc14.explicit.txt` は配布物のURL・MD5を固定した環境仕様。`linux-gcc14.packages.json` にはSHA-256も保存している。製品の数値依存とは別の、ビルド検証用ツールである。

## 環境の作成

プロジェクトルートのLinuxシェルで実行する。`micromamba` は実行可能ファイルのパスに置き換えてよい。導入先はLinuxの大文字・小文字を区別するファイルシステムに置く。

```sh
uvlas_tool_base="$HOME/.local/share/uvlas/toolchains/gcc14-cmake4.1.2"
uvlas_tool_prefix="$uvlas_tool_base/env"
CONDA_PKGS_DIRS="$uvlas_tool_base/pkgs" micromamba --no-rc -r "$uvlas_tool_base/root" create --yes -p "$uvlas_tool_prefix" --file tools/toolchains/linux-gcc14.explicit.txt
```

## 構成・ビルド・試験

```sh
micromamba --no-rc -r "$uvlas_tool_base/root" run -p "$uvlas_tool_prefix" cmake --preset linux-debug -DCMAKE_CXX_COMPILER="$uvlas_tool_prefix/bin/x86_64-conda-linux-gnu-c++"
micromamba --no-rc -r "$uvlas_tool_base/root" run -p "$uvlas_tool_prefix" cmake --build --preset linux-debug
micromamba --no-rc -r "$uvlas_tool_base/root" run -p "$uvlas_tool_prefix" ctest --preset linux-debug
micromamba --no-rc -r "$uvlas_tool_base/root" run -p "$uvlas_tool_prefix" cmake --preset linux-release -DCMAKE_CXX_COMPILER="$uvlas_tool_prefix/bin/x86_64-conda-linux-gnu-c++"
micromamba --no-rc -r "$uvlas_tool_base/root" run -p "$uvlas_tool_prefix" cmake --build --preset linux-release
micromamba --no-rc -r "$uvlas_tool_base/root" run -p "$uvlas_tool_prefix" ctest --preset linux-release
```

この環境を有効化したプロセス内でのみ実行する。生成物の実行に同環境のライブラリが必要になる場合があるため、このまま独立配布用バイナリとは扱わない。製品の配布互換性は別途検証する。

Windowsで検証したCMake 4.1.2は、プロジェクト内のPython仮想環境へ公式PyPIから導入した。システムのCMakeは変更していない。

```powershell
python -m venv build/tooling-venv
.\build\tooling-venv\Scripts\python.exe -m pip install --index-url https://pypi.org/simple cmake==4.1.2
.\build\tooling-venv\Scripts\cmake.exe --preset windows-msvc
.\build\tooling-venv\Scripts\cmake.exe --build --preset windows-debug
.\build\tooling-venv\Scripts\ctest.exe --preset windows-debug
.\build\tooling-venv\Scripts\cmake.exe --build --preset windows-release
.\build\tooling-venv\Scripts\ctest.exe --preset windows-release
```
