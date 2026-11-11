"""Report successful numerical checks and the failed adoption criterion together."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    source=F/'q9-distorted-plate-result.json';d=json.loads(source.read_text());a=json.loads((F/'q9-distorted-plate-artifact-audit.json').read_text())
    assert a['all_criteria_pass'] and a['result_sha256']==digest(source)
    lines=['# Q9の歪み格子・非一様薄肉曲げ：実施結果','','2026-10-31。**薄肉一様性の事前基準が未達で、全体判定は不合格。補正なし完全積分Q9の最終採用は保留する。** 最細の精度や最後の観測次数が成立しても、別の必須条件の不合格を打ち消さない。','','[事前仕様](q9-distorted-plate-experiment.md) の8ケースを実施した。[主結果](q9-distorted-plate-result.json) と [観測・成果物監査](q9-distorted-plate-artifact-audit.json) を保存。従来の矩形格子の結果は上書きしていない。','','## 格子収束','','表は9点組立・90桁保持解から投影した物理変位で評価した観測量。保持解の残差保証と投影状態の残差を区別する。誤差は独立な厚さ二次の面内連続参照との差。','','| 厚さ | 分割 | 中央たわみ誤差 | L2誤差 | コンプライアンス誤差 | 求積・求解等の必須条件 |','|---:|---:|---:|---:|---:|---|']
    for c in d['cases']:
        o=c['observables'];lines.append(f"| {c['h']} | {c['n']}×{c['n']} | {100*o['center_relative']:.6f}% | {100*o['l2_relative']:.6f}% | {100*o['compliance_relative']:.6f}% | {'成立' if c['eligible'] else '未達'} |")
    lines += ['','最細16×16は両厚さの3指標全て1%以内。8→16のL2観測次数は次表のとおり。全列を保持しており、この2点だけから一般の漸近性や尾部を保証しない。','','| 厚さ | 最後のL2次数 | 精度1%以内 | 次数1.5以上 |','|---:|---:|---|---|']
    for c in d['thickness_checks']:lines.append(f"| {c['h']} | {c['last_l2_order']:.6f} | {'成立' if c['accuracy_pass'] else '未達'} | {'成立' if c['order_pass'] else '未達'} |")
    mid=next(c for c in d['cases'] if c['h']==.01 and c['n']==8);thin=next(c for c in d['cases'] if c['h']==.001 and c['n']==8)
    lines += ['',f"一方、8×8での厚さ0.001のL2誤差は{thin['observables']['l2_relative']:.9g}、厚さ0.01は{mid['observables']['l2_relative']:.9g}。比は{a['thin_to_mid_l2_ratio']:.6f}で、基準の2倍＋1e−8を超えた。この歪み族で薄肉に一様な収束を支持するという採用条件は成立しない。",'','## 数値処理の確認','',f"面内7／9点の共通尺度の剛性差は最大{max(c['stiffness_7_vs_9'] for c in d['cases']):.6e}、荷重差は最大{max(c['load_7_vs_9'] for c in d['cases']):.6e}、観測量とr/d/c別応答の求積感度は最大{max(v for c in d['cases'] for v in c['quadrature_response'].values()):.6e}で、事前条件内。70／90桁の観測用物理変位は今回全ケースで同じ倍精度値になった。これは無限精度への誤差0を意味しない。",'', '保持解の残差・2種の後退誤差、対称性、正対角、全域幾何、X²補間、矩形への回帰、反転幾何の拒否、誤ったヤコビアンの負の対照が成立した。', '',f"監査は保存場から別のQ9評価・11点積分・物理点逆写像・4変数参照を再計算し、最大L2差{a['maximum_l2_difference']:.6e}で判定が一致した。監査成立と主試験不合格を両方保存する。疎行列と保持解自体は保存していないため、保持残差は記録とソースの監査で、残差全量の独立再計算ではない。",'',f"実行は{d['seconds']:.2f}秒・ピーク{d['peak_bytes']:,} bytes、監査は{a['seconds']:.2f}秒・ピーク{a['peak_bytes']:,} bytes。NPZ合計は{sum((ROOT/c['artifact']).stat().st_size for c in d['cases']):,} bytes。いずれも事前の資源条件内。",'','## 設計判断と次の作業','','数値処理の必須条件が成立した範囲で、歪みと薄肉に関する採用上の不足を追加した。係数生成・表現の全誤差を上界評価したわけではない。単一の比だけからせん断ロッキングが唯一の原因と断定しない。','','[双一次写像の無せん断適合性](../q9-distorted-shear-compatibility.md) で、Q9のたわみに必要な無せん断方向場が同じQ9へ含まれない解析的反例を整理した。この別の反例は、今回の平衡誤差の全量を説明したものではない。次は同じ物理問題のひずみ・エネルギー成分を分離し、対策する成分と補間／投影／内部変数の式を確定して候補比較へ進む。','','一般曲面、積層、梁、偏心質量、接合、支持・減衰、出力回復の採用条件は残る。J1・製品本体の開始・初回実務適用の受入は解除しない。','']
    (F/'q9-distorted-plate-results.md').write_text('\n'.join(lines),encoding='utf-8')
if __name__=='__main__':main()
