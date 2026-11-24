# Nastran の点支持と支持集合を読む

認証済みの Nastran デッキから、SPC・SPC1・SPCADD を型付きで読み取る。結果は元デッキの宣言元と SHA-256 を記録し、各欄の位置と処理に必要な値を所有する。元の座標カードは読み取るだけなので、その後もシェルなどの読取りへ渡せる。 [IO4-03](../../specs/foundation/io-execution-adoption-specification.md#io4-03-nastran)

SPC は1カードの二つの GRID 指定を読み、成分番号1〜6の集合と、入力にある場合は指定運動の正確な十進値を保持する。空欄は元のまま保持し、規格上の運動値の既定は零である。SPC1 の列挙では各 GRID の宣言を照合する。THRU では端点の値と欄を保持し、範囲内に宣言されていない番号があっても入力エラーにしない。 [MSC Nastran 2022.1 Quick Reference Guide、pp.3116–3119](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)

SPCADD の参照先は SPC または SPC1 の直接集合に限り、重複部員、未解決参照、SPCADD への入れ子を拒否する。直接集合と組合せが同じ SID を持つことは許す。両方を元欄付きで保持し、選択時に SPCADD を優先する規則は後段が適用する。Case Control でどの集合を選ぶかと、外部成分を物理的な接続・支持へ移す処理は別の段階で行う。 [MSC Nastran 2022.1 Quick Reference Guide、p.3120](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)

Windows の MSVC 19.51 Release と GCC 14.2 Release で、座標・シェル・支持を同じスナップショットから読む対象試験が各1件通過した。重複の際は二つ目の欄を指すこと、8文字を超える通常の free 整数を拒否すること、作業・確保が1不足する入力と取消で元デッキを保持することを照合した。clang-tidy 21.1.1 は読取りと試験を検査した。隔離した GCC コピーで入れ子の拒否を外すと回帰が失敗し、復元後は通過した。
