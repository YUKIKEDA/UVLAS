# MAT8 の面内定数と PCOMP の宣言層を読む

SHA-256 照合済みの Nastran デッキから、MAT8 の定数と PCOMP の積層を型付きで取り出す。結果はシェルの PID から PCOMP、各層の MID から MAT1 または MAT8 の宣言と元の欄をたどれる。[IO4-03](../../specs/foundation/io-execution-adoption-specification.md#io4-03-nastran)

## MAT8 の定数を入力と同じ精度で保つ

E1、E2、NU12 を読み、明示した弾性率は正、G1Z と G2Z は非負とする。G12、G1Z、G2Z、RHO の空欄と明示値を区別する。熱膨張係数、参照温度、強度欄、減衰、F12、STRN も保持する。STRN の明示値は十進値として厳密に 1.0 と照合する。Hashin の追加行は元欄を示して未対応と返す。面外の材料定数やせん断剛性の補完、単位の確定は後段で行う。[MSC Nastran 2022.1 Quick Reference Guide、pp.2185–2188](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)

## PCOMP は宣言された層順と積層モードを持つ

層を底面から宣言順に保持する。MID、厚さ、角度、SOUT の4欄がすべて空の組は層を作らない。最初の層には MID と厚さが必要で、後続の MID と厚さが空欄なら、それぞれを最後に明示した層の番号を保持する。厚さの値や材料欄を埋めずに出所をたどれる。角度は要素材料軸からの度の指定として読み、SOUT の YES、NO、空欄を区別する。

LAM は空欄の全積層、SYM の対称半積層、MEM の膜のみ、BEND の曲げのみ、SMEAR の均質化指定、SMCORE のコア付き指定を区別する。SYM の展開や SMCORE の厚さ分配は後段で行う。Z0 は基準面から層の底面までの距離として、要素の ZOFFS とは別の欄で保持する。非構造質量、接着層の許容せん断応力、参照温度、減衰も元値で保持する。FT の語形は採用版の7種類に限る。いずれの指定にも SB（接着層の許容せん断応力）が必要である。未知の LAM、FT、SOUT は欄位置を示して拒否する。総厚、既定 Z0、破壊指数の評価、3D 材料の適格性は後段で扱う。[MSC Nastran 2022.1 Quick Reference Guide、pp.2764–2768](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf)

## 検証

同じスナップショットを字句、座標、シェル、材料、積層の順に渡し、元バイトと SHA-256 を照合した。small・large・free の PCOMP で層順、空の組、継承元、選択肢を確認した。MAT8 の十進値と欄位置、未知オプション、未解決の材料、最初の層の欠落、メモリと作業の境界、取消を試験した。

Windows の MSVC 19.51.36260.0 Debug と Linux の GCC 14.2.0 Debug（ASan・UBSan）で、Nastran の対象6件がそれぞれ通過した。積層カードは各5ケース・413アサーション。Linux の clang-tidy 21.1.1 で本体の通常構成・確保計測構成と試験を検査した。
