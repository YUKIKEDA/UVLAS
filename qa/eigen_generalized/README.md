# 対称一般化スペクトルの基準結果

製品の一般化抽出で次の固定行列を解く。

$$
M=\operatorname{diag}(1,4,16),\qquad K=\operatorname{diag}(-\delta,0,16\delta),\qquad \delta=2^{-60}.
$$

解析根は $$ -\delta,0,\delta $$、質量正規化basisは $$ \operatorname{diag}(1,1/2,1/4) $$ である。尺度は $$ D=I,s_\lambda=m_0=1 $$ と固定する。

全3根、basisの対角成分、元残差・M正規化・Gram誤差の保守上界を16進浮動小数点で保存する。微小非零値を零へ丸めず、負・零・正を区別する。独立試験では非対角Mと非2冪尺度も使い、有理数で原問題を照合する。

基準は`tools/qa/results.py accept --case eigen-generalized --run <executable>`で実測から登録する。別OSの値を推測して登録しない。CIだけにあるplatformは実ログを`--output`へ渡す。この小系の成功は物理modeの被覆や実務性能の受入を意味しない。
