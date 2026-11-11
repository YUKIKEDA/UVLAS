"""Freeze MR-R1 closure and shared-space rejection without replacing old records."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def read(n):return json.loads((F/n).read_text())
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    reference=read('plate-modal-exact-reference-result.json');evaluation=read('plate-modal-reference-reevaluation.json');space=read('q9-shear-gradient-compatibility-result.json');audit=read('q9-shear-gradient-compatibility-audit.json')
    assert all(x['all_criteria_pass'] for x in [reference,evaluation,space,audit])
    assert evaluation['reference_sha256']==digest(F/'plate-modal-exact-reference-result.json') and audit['parent_sha256']==digest(F/'q9-shear-gradient-compatibility-result.json')
    lines=['# 4変数参照の精度確定と共有せん断空間の審査結果','','2026-10-31。**参照の数値精度問題を切り分け、MR-R1による独立観測監査を完了した。一方、P-S36は今回の歪み格子で共有後の無せん断たわみ空間を失うことを確認し、製品採用を見送る判断を維持する。**','','## MR-R1の精度と観測監査','','[参照精度契約](../plate-modal-reference-accuracy-contract.md) で入力h・π・荷重・実際の弾性行列のbinary64値を固定した。厚さ多項式の有理数積分・正確な4×4消去と、別の閉形式によるDecimal80／100桁LDLを照合した。行列閉形式と正確な残差0が成立し、最近接binary64への丸めは0.5ULP以内、別経路からのbit一致を確認した。','','| 厚さ | 旧中央振幅の相対差 | Decimal100と有理数解の最大相対差 | π入力感度（中央振幅） |','|---:|---:|---:|---:|']
    for c in reference['cases']:lines.append(f"| {c['h']} | {c['old_center_relative_difference']:.6e} | {c['decimal_paths']['100']['maximum_relative_difference']} | {c['pi_input_center_sensitivity']:.6e} |")
    lines += ['', 'πの構造的な入力感度は別のMachin級数評価で診断した。固定binary64入力に対する代数精度と、真のπ・連続体・厚さモデルへの差を区別する。', '', '[保存場の再評価](plate-modal-reference-reevaluation.json) は元Q9とP-S36の各8ケース、7／9点組立×70／90桁の保存場を対象にした。主経路は有理数解、別観測経路はDecimal100の解から同じ参照を丸め、別形状関数・物理中心探索・11点積分で再計算した。中心／仕事／L2／最後の次数と参照丸めの伝播予算が全て成立。有限要素の再求解は行っていない。','','| 方式 | MR-R1独立観測監査 | n=8の薄肉／中厚L2比 | 薄肉一様性・全体受入 |','|---|---|---:|---|']
    for d in evaluation['datasets']:lines.append(f"| {d['name']} | {'成立' if d['independent_observation_audit_pass'] else '未達'} | {d['thin_to_mid_l2_ratio']:.6f} | {'成立' if d['experiment_acceptance_pass'] else '未達'} |")
    lines += ['', '旧5／7点の倍精度参照・旧監査false・元の不合格結果を保持した。新しい参照版での監査成立であり、旧結果を後付けで書き換えたものではない。両方式の薄肉一様性が未達という判断は変わらない。', '', '## 局所から共有自由度への適合性','','[適合性契約](../q9-shear-gradient-compatibility-contract.md) に従い、厚さ0次の12せん断拘束を取り出した。局所では方向場rank12で、局所Q9たわみに対応する方向場が存在した。しかし同じ物理節点と支持条件を共有すると、次の結果となった。','','| 分割 | 矩形の許容たわみ次元 | 歪み0.2の許容たわみ次元 |','|---:|---:|---:|']
    for n in [1,2,4,8]:
        a=next(c for c in space['cases'] if c['n']==n and c['amplitude']==0);b=next(c for c in space['cases'] if c['n']==n and c['amplitude']==.2)
        lines.append(f"| {n}×{n} | {a['ranks']['1e-10']['compatible_w_dimension']} | {b['ranks']['1e-10']['compatible_w_dimension']} |")
    lines += ['', 'n=1は歪みの内部頂点がなく同一格子。SVDの3閾値、左零空間の欠損作用素、ピボットQRのrankが一致した。[別の2D積分監査](q9-shear-gradient-compatibility-audit.json) でも作用素と次元が一致し、D,Wの最大相対差は約1.2e−15だった。', '', '[理論的な解釈と条件付き下界](../q9-shared-shear-space-disposition.md) では、厚さ二次係数を含む全せん断場について、二次係数を自由にしても厚さ0次の欠損を完全には消せないことを示した。固定格子での許容空間の喪失を、局所ゲートの良好な結果で相殺しない。これは今回の有限作用素に対する診断で、一般曲面の連続積分や全誤差の原因を厳密に証明した意味ではない。', '', '## 現在の決定','','- この限定4変数参照の数値入力・精度・丸め・観測への伝播はMR-R1として閉じた。NU-02全体や実務実験の参照採用を閉じたことにはしない。','- P-S36の現構成は、今回の共有・境界条件付き歪み格子で必要な無せん断空間を保持しない。製品対策として採用せず、このままの追加細分化を次作業にはしない。','- 次の候補は辺の物理トレースを保存する共変ひずみ補間を審査する。還元と勾配の可換性、隣接要素との対応、零モード、材料と曲面の整合を先に定め、共有後の適合性ゲートを大域比較の前に通す。','- 構造モデル全体、J1、主要製品I/F、実装開始条件は未完了のまま。','']
    path=F/'modal-reference-and-shear-space-results.md';path.write_text('\n'.join(lines),encoding='utf-8')
    names=['plate-modal-exact-reference-result.json','plate-modal-reference-reevaluation.json','q9-shear-gradient-compatibility-result.json','q9-shear-gradient-compatibility-audit.json','modal-reference-and-shear-space-results.md']
    result=dict(source_sha256=digest(__file__),evidence_sha256={n:digest(F/n) for n in names},modal_reference_version='MR-R1',modal_reference_accuracy_closed=True,
        independent_observation_audits_pass=True,P_S36_adopted=False,shared_shear_space_diagnostic_completed=True,product_shell_adopted=False)
    (F/'modal-reference-and-shear-space-disposition.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
