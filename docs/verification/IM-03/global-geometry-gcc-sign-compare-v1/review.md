# GCCの符号比較診断への修正

PR #153のLinux Debug・Release・Sanitizerは、alias_kinematics.cppの同じ符号比較でbuildに失敗した。個別job 134885864495の実ログを確認した。uint16_tのfirst_limbは整数昇格でintとなり、256-first_limbもintになる。一方、span::size()はsize_tで、GCC14の-Werror=sign-compareが拒否した。

上限をstd::size_t{256}として減算し、両辺を同じ非負寸法の型にした。直前のfirst_limb>=256による短絡検査は保持した。受理する入力、資源予算、保持値と数式は変わらない。新製品4TUと3内部headerのsize比較・小整数型・減算を検索し、同類はこの1箇所だった。追加依頼のdetail()再取得型も新試験2TUへ限定検索し、該当しなかった。

専用Debugで大域幾何normal/resourcesの2targetだけをbuildし、CTestは2/2成功した。Clang21.1.1で全新4製品TUに-Wsign-conversion/-Wconversion/-Wshadow/-Wsign-compare/-Werror=sign-compareとbugprone-narrowing-conversionsを適用し、exit0だった。出力の76 warningsは既存外部依存由来で、製品の診断はなかった。ClangとGCCの警告判定は同一とは限らず、このWindows検査をLinux成功の証拠にはしない。

基点は82de6ad。既存の固定source・変異・測定記録は書き換えていない。今回の差分は上限式の型指定だけなので、新しい数値方式試験や隔離変異、無関係suite、Release全体の再実行は不要と判断した。Linux3構成はpush後のCIへ委ね、待機しない。manifestは訂正後sourceと検査出力のSHA256を記録する。
