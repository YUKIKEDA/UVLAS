# SI変換の隔離変異検証

固定製品source `5eaaeed16771ed630463348343f1433cad696efc` で、controlと6件の意味変異を別々にcompile/linkして実行した。controlはnormal 17 cases / 39,034 assertions、interval 4 cases / 55 assertionsが成功。6変異はいずれもcompile/link成功後にCatchのassertionで失敗した（exit 42）。compile失敗を検出に数えていない。

| 変異 | 検出した独立条件 |
|---|---|
| length-factor-inverse | mmの倍率を逆にし、独立有理値37/1000の包含が失敗 |
| degree-pi-factor-replaced | 単位のπ冪を1にし、Support.elasticの独立区間包含が失敗 |
| two-word-low-dropped | exact fallbackの公開TwoWord lowを0にし、疎な低位語を持つ有理値との一致が失敗 |
| foreign-catalog-checks-disabled | requestの早期owner検査と後段の完全key比較を両方抑止し、他catalog由来keyの拒否が失敗 |
| work-meter-disabled | Work::stepの使用量加算を除き、独立の予定量5に対する実行量0を検出 |
| division-directed-rounding | 除算の外向き丸め方向を反転し、独立有理数1/3の包含が失敗 |

normal4TUの既存test objects、IO/Structure/CoreのDebug libを、実行開始時に私有supplyへ複製した。供給時点は同sourceのclean Debug build直後。全6製品TUを私有sourceで再compileし、archiveより先に明示objectをlinkした。inline Work/丸めの変異は私有headerを含むinterval test TUを再compileした。供給hashは`supply-sha256.json`、全sourceのgit blobは`summary.json`にある。

初回controlはstructure内部headerのコピー不足でcompile停止した。コピー範囲を修正し、git archiveで同SHAのinclude/src/io/src/structure/tests/ioを展開してから成功controlを実行した。作業treeで後続のOpcode改名が始まったため、現作業treeは読まず、この固定展開と開始時のsupplyを使った。πのanchorが重複、次いで改行差で未一致となったproducer停止も、変異検出に数えていない。成功済みcontrol/変異は同じ隔離領域で保持し、Resumeで残りだけを実行した。最終producerのexitは0。

foreign-keyの早期owns_key検査だけを除く予備変異はCatch 0だった。後段のrow.key==entry.keyが同じ世代を検査するためである。その成功ログ・変異sourceも別名foreign-catalog-generationで保持し、6件の検出成功には混ぜていない。製品の冗長検査を残したまま、両guardを無効にする意味変異で負対照の感度を確認した。

各Catch実行のseed、最初の不一致、compile/link出力は`logs-and-mutants.zip`にある。ログは個人パスを置換し、行末空白を除いたLF表示記録であり、生stdoutのbyte同一を主張しない。検出試験は固定有理値・固定入力で、乱数の標本抽出は使わない。再現時はログのseedをCatchの--rng-seedへ渡せる。

再実行には新しいscratchを用意し、run-isolated.ps1をその直下へコピーする。同SHAをgit archiveでscratch/issue148-fixed-sourceへ展開する。対応するclean Debug buildの3 libと4 normal test objectsをscratch/issue148-mutations/supplyへコピーし、供給SHA256を照合する。元checkoutのout/issue148-build/compile_commands.jsonを保ち、MSVC環境で -Root <checkout> -Output <scratch/new-result> -SourceCommit 5eaaeed... を指定する。出力の上書きは拒否し、Resumeは中断された同じsourceの完了caseだけを再利用する。Windows MSVC19.51と固定vcpkg prefixに対する検証producerであり、Linux実行は主張しない。

これは固定sourceにおける試験感度の記録であり、他の資源故障点の成功、後続の製品修正、全体ソルバーの認定を代替しない。製品本体・既存証拠は変更していない。
