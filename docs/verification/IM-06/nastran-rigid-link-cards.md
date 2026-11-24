# Nastran の RBE2 の節点と熱膨張欄を読む

要素 ID の重複を検査したデッキから RBE2 を読む。独立節点 GN と宣言順の従属節点 GM を GRID と照合し、従属自由度 CM をビット集合にする。元デッキの SHA-256、宣言元、各欄の位置を保持する。 [IO4-03](../../specs/foundation/io-execution-adoption-specification.md#io4-03-nastran)

GM の列に実数表記が現れたら熱膨張係数 ALPHA として読み、次の実数欄を基準温度 TREF として読む。数値の十進表現と省略の区別を保持する。ALPHA 欄だけを空けて数値の TREF を後置する形式は、この読取りの対象外とする。カード間の従属自由度や拘束ループの検査、拘束式の生成、熱膨張の適用はモデル準備に残る。 [MSC Nastran 2022.1 Quick Reference Guide、pp.2976–2977](https://help-be.hexagonmi.com/bundle/MSC_Nastran_2022.1_Quick_Reference_Guide/raw/resource/enus/MSC_Nastran_2022.1_Quick_Reference_Guide.pdf#page=3012)

Windows の MSVC 19.51 Release と Linux の GCC 14.2 Release で、自由形式と固定幅継続、節点順、自由度、ALPHA と TREF の元位置を確認した。未定義節点、節点の重複、不正な自由度、余分な熱膨張欄も試した。clang-tidy 21.1.1 は読取りと試験を検査した。GCC の分離コピーで TREF を ALPHA の欄から読むように変えると対象試験が失敗し、復元後は通過した。
