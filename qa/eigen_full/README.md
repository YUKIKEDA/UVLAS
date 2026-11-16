# 完全スペクトルと派生結果の基準

製品の完全法を通し、全固有値・原Schur認定・単一pair・左右clusterの保存値をQA-03で照合する。入力は $$ A=\operatorname{diag}(-2^{-60},2,3) $$ 、成分尺度と固有値尺度は1で固定する。負の微小根は零へ置き換えない。

出力は3固有値、原Schur残差・Gram誤差、固有値2に対応するpairの正規化誤差・右vector第2成分、同じ1次元clusterの左原残差の計8値。十進に丸めず16進で出力する。残差と正規化の保存値は製品の保守上界であり、backendが返した成功状態だけを保存するものではない。

`tools/qa/results.py accept --case eigen-full --run <executable>` で、cleanな製品sourceからプラットフォーム別基準を作る。別OSの値をコピーして登録しない。`qa.eigen_full` は登録値との比較を行う。

一般の非正規行列、欠損固有basis、非連続cluster、非2冪尺度、PG24の原作用は専用数値試験で照合する。この小さい基準caseを物理mode被覆・実務受入・製品規模の性能の証拠とはしない。
