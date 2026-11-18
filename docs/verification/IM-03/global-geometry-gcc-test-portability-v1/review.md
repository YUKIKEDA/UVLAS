# 大域幾何試験のGCC移植性修正

c94025aのLinux3構成で、試験側に同じ4診断が出た。3つのpair range-loopをconst参照にし、vector::reserveへ渡すworker系列と対応するloop添字をsize_tへ揃えた。入力値、全組合せ、assertion、許容差は変わらない。製品本体は変更していない。

対象試験2TUのrange-loopを検索し、残る値渡しはscalarだけと確認した。Windows専用Debugで2targetをbuildし、CTest 2/2成功。Clang21.1.1の通常clang-tidyに-Wsign-conversion/-Wconversion/-Wrange-loop-construct/-Werror=range-loop-construct/-Wsign-compareを加えて2TUを検査し、exit0だった。外部header等の抑制された診断数を、製品診断やLinux成功と取り違えない。formatも成功した。

今回の試験側訂正は同じ入力・呼出しを維持するため、Release全体、隔離変異、無関係suiteは再実行しない。旧証拠は保持し、この追補だけを追加する。LinuxとWindowsのCIはpush後に動作し、ここでは完了待ちをしない。Windowsのinput_document.cpp診断は親担当が別に調べており、この修正では触れていない。
