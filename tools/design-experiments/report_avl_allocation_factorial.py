"""Render reviewed factorial evidence without treating sensitivity as error."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def read(n):return json.loads((F/n).read_text())
def main():
    runs=read('avl-allocation-factorial-runs.json');data=read('avl-allocation-factorial-result.json');audit=read('avl-allocation-factorial-artifact-audit.json')
    assert runs['completed'] and data['methods_criteria_pass'] and audit['all_criteria_pass']
    lines=['# 内側・翼端の独立配分対照：実施結果','','2026-10-31。[事前仕様](avl-allocation-factorial-experiment.md) の新規4計算を実施した。既存A（一様）・D（tip-rich）の元結果を維持し、B（内側のみ変更）・C（翼端のみ変更）を追加した。以下の差は有限格子の配分感度であり、連続解の誤差や精度改善率ではない。','','## 実行と検査','','| 翼 | 変更 | 未知数 | 秒 | ピークRSS（bytes） |','|---|---|---:|---:|---:|']
    for r in runs['runs']:lines.append(f"| {r['case']} | {r['allocation']} | {r['native_unknowns']} | {r['seconds']:.2f} | {r['peak_resident_bytes']} |")
    lines += ['',f"4計算の実行段階は{runs['seconds']:.2f}秒、最大RSSは{max(r['peak_resident_bytes'] for r in runs['runs']):,} bytes。観測前後の8出力一致、after ZIPの無損失保存、入力・実行体・元結果のハッシュを確認した。後処理は{data['seconds']:.2f}秒、ピーク{data['peak_bytes']:,} bytes、今回の出力と派生場の合計{data['combined_output_bytes']:,} bytes。",'', '8場の実区間別帯数、残差1e−11以内、全域荷重と20帯の保存1e−10以内、配分差の再構成1e−12以内、地域別二乗ノルム分割、事前資源上限を確認した。独立したDecimal60の帯データ再構成でも効果ベクトルとノルムを照合した。これは独立な数値再計算であり、人による独立審査や全係数の独立検証を意味しない。','','## 内側・翼端・交互作用','','A/B/C/Dは事前仕様の4頂点。内側効果B−A、翼端効果C−A、交互作用D−B−C＋Aの和を、既存のD−Aと照合した。表は各効果ノルム／D−Aノルム。ノルムには交差項があり、足して1になる寄与率ではない。','','| 翼 | 観測量 | 内側／合成 | 翼端／合成 | 交互作用／合成 | 内側と合成の方向余弦 | 翼端と合成の方向余弦 |','|---|---|---:|---:|---:|---:|---:|']
    labels={'tau0':'τ帯平均','tau_p1':'τ帯P1','k0':'K帯平均','j_line':'有限辺モーメント','j_point':'元点モーメント'}
    for response in data['responses']:
        for channel,item in response['channels'].items():
            e=item['effects'];den=e['combined']['norm'];ratio=lambda n:e[n]['norm']/den
            lines.append(f"| {response['case']} | {labels[channel]} | {ratio('inner'):.6f} | {ratio('tip'):.6f} | {ratio('interaction'):.6f} | {e['inner']['cosine_to_combined']:.6f} | {e['tip']['cosine_to_combined']:.6f} |")
    lines += ['','## 全域荷重傾斜の変化','','各列は対応する有限格子差の符号付き値。分布のノルムと区別し、総量の相殺を分布精度の根拠としない。','','| 翼 | 量 | 内側 | 翼端 | 交互作用 | 合成 |','|---|---|---:|---:|---:|---:|']
    for response in data['responses']:
        for channel,item in response['scalars'].items():
            values=[item['effects'][n]['signed']['total'] for n in ['inner','tip','interaction','combined']]
            lines.append('| '+response['case']+' | '+channel+' | '+' | '.join(f'{v:.9e}' for v in values)+' |')
    lines += ['','## 判定','','この4頂点では内側と翼端の変更、および交互作用を区別できた。帯ごとの差、地域別二乗ノルム、各効果と合成の方向余弦、交差項、全域荷重を [数値結果](avl-allocation-factorial-result.json) に保存した。[成果物監査](avl-allocation-factorial-artifact-audit.json) は原本・ZIP・資源・実配分・Decimal再計算を確認した。','','非相似格子なのでRichardson次数・連続解への尾部・精度改善率は求めない。M12の診断からM24の誤差予算を補完しない。今回の検査成立だけでAVLを製品精度の受入参照に採用しない。現在の採否と後続設計への影響は [AVL診断の終了判定](avl-reference-diagnostic-disposition.md) にまとめる。','']
    (F/'avl-allocation-factorial-results.md').write_text('\n'.join(lines),encoding='utf-8')
if __name__=='__main__':main()
