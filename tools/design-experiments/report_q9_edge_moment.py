"""Record constructive feasibility and its remaining adoption boundaries."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[2];F=ROOT/'docs/verification'
def read(n):return json.loads((F/n).read_text())
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    data=read('q9-edge-moment-feasibility-result.json');audit=read('q9-edge-moment-feasibility-audit.json');assert data['all_criteria_pass'] and audit['all_criteria_pass'];assert audit['parent_sha256']==digest(F/'q9-edge-moment-feasibility-result.json')
    lines=['# 幾何依存の辺再現と仕事補正の可解性結果','','2026-11-02。**G-EM1は今回の40配置で一定物理ひずみ・一定応力の節点仕事を保持し、10共有辺対でトレースを保持した。自由要素rank75と平坦T-S36の継承も成立した。** 一般曲面の共有後の無せん断空間・非線形接線・材料・初期ひずみの製品採用は未完了。','','## 構成した作用素','','[候補の定式化](../q9-edge-moment-candidate.md) は次の3段からなる。','','1. 共有辺の接線・厚さ線だけから、一定物理テンソルの観測可能な成分を抽出する。未観測成分を対角正則化で作らず、Mandel内積の最小ノルム解を使う。','2. 一定場の体積再現欠損を補い、隣接要素で異なる体積平均の影響を辺上で打ち消す。辺の補正は共有情報だけに依存させる。','3. 辺で0、体積平均1のbubbleで残る体積平均差を除く。一定場では補正0となり、共有トレースを変えずに一定応力の仕事を保持する。','','前回のG-COV/C-COVをそのまま採用したものではなく、値と仕事のどちらかだけを修正したものでもない。辺の幾何依存再現を先に作ることで、内部bubbleだけでは直せなかった辺上欠損を処理した。','','## 事前ゲートと結果','','[事前仕様](q9-edge-moment-feasibility-experiment.md) のRECT/SKEW/TRAP/SADDLE/NORMAL-C、厚さ0.1／0.01、剛体回転・同じ物理厚さの座標再表現を組み合わせた40配置を実施した。荷重平衡は解いていない。値差は一定物理Mandel6基底、仕事差は全81仮想変位に対する6一定応力の作用素差で、変位応答の誤差率ではない。','','| 検査 | 最大値または結果 | 事前基準 |','|---|---:|---:|',
        f"| 一定物理ひずみの値差（内部・辺・追加点） | {max(c['value_error'] for c in data['cases']):.6e} | 1e−10 |",
        f"| 11点再積分での一定応力の仕事差 | {max(c['work_error_11'] for c in data['cases']):.6e} | 1e−10 |",
        f"| 10共有辺対のトレース差 | {max(c['trace_error'] for c in data['shared_pairs']):.6e} | 1e−10 |",
        f"| 7／9点で構成した作用素の場の差 | {max(c['field_7_9_difference'] for c in data['cases']):.6e} | 1e−8 |",
        f"| 平坦T-S36との尺度化剛性差 | {max(c['flat_difference'] for c in data['cases'] if c['flat_difference'] is not None):.6e} | 1e−10 |",
        f"| 剛体6モードの作用 | {max(c['rigid_error'] for c in data['cases']):.6e} | 1e−11 |",
        '| 自由要素rank | 全40配置・3閾値で75 | 75 |']
    lines += ['', '辺の観測可能rankは、平坦辺で3、今回の曲面辺で6となり、3つの閾値で安定した。保持空間の条件数は最大約89.44で上限1e8以内。捨てたモードの辺上の値への影響も事前条件内だった。回転・厚さ再表現の値・仕事・剛性・幾何、アフィン場、正半定値と対称性のゲートも全て成立。', '', '仕事の補正は、構成に使った求積上では代数的に成立する。この0残差だけを根拠とせず、9点構成の作用素を11点で再積分して同じ基準を確認した。平坦剛性の継承は今回の線形作用素の一致であり、一般曲面の薄肉一様性の証拠へ移さない。','','## 別実装による監査','','[監査結果](q9-edge-moment-feasibility-audit.json) は回転・再表現前の全10基準配置と10共有辺対で成立した。一次実装の3×3テンソル処理に対し、別の形状微分と6×6 Mandel変換を使って作用素を作り直した。辺係数は手組みSVD解に対して別の最小二乗ルーチンで再計算し、保存された係数・体積平均・モーメントを照合した。','','11点の仕事だけでなく、元の一定応力の節点力は6境界面のσn積分とも照合した。独立に再構成した剛性の尺度化差は最大約9.19e−16、11点での仕事差は約5.42e−15以内。人による独立審査や全曲面の収束証明とは区別する。', '', f"主計算は{data['seconds']:.2f}秒・ピーク{data['peak_bytes']:,} bytes、監査は{audit['seconds']:.2f}秒・ピーク{audit['peak_bytes']:,} bytesで事前条件内。",'', '## 閉じた範囲と次の決定','','値・仕事・辺の3条件を同時に満たす幾何依存の辺還元／内部拡張が、今回の線形化した有限空間で可解であることを確認した。局所の安定性ゲートも成立したため、G-EM1を次の審査へ残す。旧G-COV/C-COVとP-S36の不採用、平坦T-S36の限定通過は維持する。','','次は次の事項を閉じる。','','- 試作は線形Bと一定6場への写像に限定している。一般のひずみ場を受けるRと双対R*について、[場と双対の処理契約](../q9-edge-moment-field-interface.md) に評価点・辺データ・作用順・所有権・単位・キャッシュ・失敗保証を定義した。次はこの一般場処理を実装して双対性を検査し、非線形Greenひずみへ接続する。','- 適合初期ひずみの自然状態、全構成行列・材料軸、エネルギーからの内力と幾何接線を照合する。','- 曲面で要素を共有した後の許容空間と、曲げ・膜・厚さ成分の誤差を審査する。自由1要素のrank75を、共有後の安定性や薄肉一様性の代用にしない。','','これらを閉じる前に大域荷重計算や製品APIを追加しない。主シェル最終採用・構造審査全体・J1・製品コア開始は未完了。','']
    path=F/'q9-edge-moment-results.md';path.write_text('\n'.join(lines),encoding='utf-8')
    names=['q9-edge-moment-feasibility-result.json','q9-edge-moment-feasibility-audit.json','q9-edge-moment-results.md']
    state=dict(source_sha256=digest(__file__),evidence_sha256={n:digest(F/n) for n in names},candidate='G-EM1',linear_value_work_trace_feasible=True,local_rank_gates_pass=True,flat_T_S36_inheritance_pass=True,
        independent_reconstruction_pass=True,curved_shared_shear_space_closed=False,nonlinear_material_initial_strain_closed=False,decision='retain candidate for general-field adjoint and nonlinear review',product_shell_adopted=False)
    (F/'q9-edge-moment-disposition.json').write_text(json.dumps(state,indent=2)+'\n')
if __name__=='__main__':main()
