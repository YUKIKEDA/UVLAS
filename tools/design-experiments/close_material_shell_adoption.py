"""Audit material 1.2 integration and publish adoption without changing R4/R5."""
from pathlib import Path
import hashlib,json,re,time
import numpy as np
import material_equivalent_input as inp
import check_material_shell_integration_v2 as checker
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(x,message):
    if not x:raise AssertionError(message)
def main():
    tick=time.monotonic();dest=ROOT/'docs/verification/material-shell-adoption-disposition.json';require(not dest.exists(),'immutable output')
    r5path=ROOT/'docs/verification/paired-shell-r5-audit.json';r5=json.loads(r5path.read_text());c5path=ROOT/'docs/verification/paired-shell-r5-closure.json';c5=json.loads(c5path.read_text())
    require(r5['passed'] and c5['main_shell_adopted'],'adopted main shell')
    sources={};artifacts={}
    for name,h in c5['immutable'].items():require(sha(ROOT/name)==h,'changed R5 '+name);sources[name]=h
    for label in ['sources','artifacts']:
        for name,h in r5[label].items():require(sha(ROOT/name)==h,'changed predecessor '+name)
        (sources if label=='sources' else artifacts).update(r5[label])
    sources[c5path.relative_to(ROOT).as_posix()]=sha(c5path)
    path=ROOT/'docs/verification/material-shell-integration-v2-result.json';r=json.loads(path.read_text());require(r['passed'] and len(r['records'])==10 and not r['failures'],'new integration gates')
    bycase={q['case']:q for q in r['records']}
    require(bycase['invalid_inputs']['count']==28 and len(set(bycase['invalid_inputs']['rejected']))==28,'negative controls')
    for q in r['records']:
        if 'metrics' in q:require(max(q['metrics'].values())<=1e-8,'material metric '+q['case'])
    require(max(bycase['quadratic_strain_velocity_section_work']['metrics'].values())<=1e-10,'section work')
    require(bycase['equivalent_density_and_physical_spatial_inertia']['inertia_error']<=1e-9,'physical inertia')
    require(bycase['nonunique_ABD_and_areal_mass']['normal_energy_difference']>1e-4 and bycase['nonunique_ABD_and_areal_mass']['second_mass_moment_difference']>1e-7,'insufficient input counterexamples')
    require(bycase['normal_and_material_basis_reversal']['tensor_error']<=1e-12 and bycase['normal_and_material_basis_reversal']['forgotten_component_transform_defect']>1e-3,'basis reversal')
    require(all(q['stress_difference']>1e-6 and q['lower']!=q['upper'] for q in bycase['unequal_layers_two_sided_recovery']['interfaces']),'two-sided recovery')
    for reportpath in [path,ROOT/'docs/verification/material-shell-integration-result.json',ROOT/'docs/verification/material-layer-kernel-result.json']:
        report=json.loads(reportpath.read_text());sources[reportpath.relative_to(ROOT).as_posix()]=sha(reportpath)
        for p,h in report['sources'].items():require(sha(ROOT/p)==h,'changed source '+p);sources[p]=h
    require(not json.loads((ROOT/'docs/verification/material-shell-integration-result.json').read_text())['passed'],'initial adapter failure retained')
    require(json.loads((ROOT/'docs/verification/material-layer-kernel-result.json').read_text())['all_checks_pass'],'legacy constitutive gates')
    artifact=ROOT/r['artifact']['artifact'];require(sha(artifact)==r['artifact']['artifact_sha256'],'snapshot');artifacts[artifact.relative_to(ROOT).as_posix()]=sha(artifact)
    a=np.load(artifact);require(a['section_K'].shape==(18,18) and a['shell_K'].shape==(513,513) and a['fullM'].shape==(513,513),'physical and section operator dimensions')
    section=np.block([[a['section_C'][i+j] for j in range(3)] for i in range(3)])
    require(np.linalg.norm(a['section_K']-section)<=1e-14*np.linalg.norm(section),'saved section matrix')
    require(np.linalg.eigvalsh(a['fullM']).min()>0 and np.linalg.norm(a['shell_K']-a['shell_K'].T)<=1e-10*np.linalg.norm(a['shell_K']),'mass and tangent')
    # Independent engineering tensor section calculation, re-evaluated against saved moments.
    req=r['material_input'];stack=inp.prepare_stack(req);rng=np.random.default_rng(20261105);q=rng.normal(size=(3,6))*np.array([.003,.04,.2])[:,None];v=rng.normal(size=(3,3))
    U,f,ke=checker.direct_section(req,.05,q,v);_,_,_,w0=inp.section_moments(stack,.05)
    expected=q.ravel()@a['section_K']@q.ravel()/2-q.ravel()@a['section_g'].ravel()+w0
    mass=np.kron(np.array([[a['section_mass'][i+j] for j in range(3)] for i in range(3)]),np.eye(3))
    repeated=dict(energy=abs(expected-U)/abs(U),force=float(np.linalg.norm(a['section_K']@q.ravel()-a['section_g'].ravel()-f.ravel())/np.linalg.norm(f)),
        mass=abs(v.ravel()@mass@v.ravel()/2-ke)/abs(ke))
    require(max(repeated.values())<=1e-10,'independent saved section work')
    # Complete the normal-reversal identity for the initial-strain vector, too.
    theta=.37;A=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1.]])
    Q=np.diag([1.,-1.,-1.]);T=inp.transform(A);Tq=inp.transform(Q);Tp=inp.transform(A@Q);m=stack[2].material
    e0error=float(np.linalg.norm(Tp@(Tq.T@m.e0)-T@m.e0)/np.linalg.norm(m.e0));require(e0error<=1e-12,'initial strain basis reversal')
    spec=ROOT/'docs/material-shell-adopted-specification.md';text=spec.read_text(encoding='utf-8');decisions=[f'MAT-D{i:02d}' for i in range(8,14)]
    require(all(d in text for d in decisions),'decisions');require(all(s in text for s in ['STR-01','1.3','ABD','MaterialUseDomain','確定状態','g_k','m_k']),'required specification boundaries')
    for p in [spec,Path(__file__),ROOT/'docs/product-requirements.md',ROOT/'docs/acceptance-criteria.md',Path(inp.__file__),Path(checker.__file__)]:sources[p.relative_to(ROOT).as_posix()]=sha(p)
    metrics=dict(integration_groups=10,negative_controls=28,section_work_error=max(repeated.values()),shell_tangent_error=bycase['three_input_kinds_same_adopted_S']['metrics']['tangent_direction'],
        shell_independent_tensor_force_error=bycase['three_input_kinds_same_adopted_S']['metrics']['independent_tensor_force'],physical_inertia_error=bycase['equivalent_density_and_physical_spatial_inertia']['inertia_error'],initial_strain_basis_reversal_error=e0error)
    disposition=dict(material_1_2_complete=True,main_shell_adopted=True,structural_review_complete=False,product_implementation_start=False,
        decisions=decisions,adopted_inputs=['isotropic 3D','orthotropic 3D nine constants','explicit full SPD 3D effective single layer or layer stack'],
        abd_only_input='reject as insufficient 3D closure; explicit completion belongs to the input model',normal_stress_condensation_in_main_shell=False,
        real_material_validation=False,material_validity_owner='work 5: data, use conditions, strain envelope and independent acceptance',
        specification=spec.relative_to(ROOT).as_posix(),metrics=metrics,next_work='1.3 beam adoption and shell interface contract')
    dest.write_text(json.dumps(disposition,indent=2),encoding='utf-8');sources[dest.relative_to(ROOT).as_posix()]=sha(dest)
    report=ROOT/'docs/verification/material-shell-adoption-results.md'
    report.write_text(f'''# 構造1.2完了：材料・積層・等価入力の採用

2026-11-05。**1.2の数値仕様・主シェル結合を採用し、次の主作業を1.3梁の採用へ移す。** [採用仕様v1.0](../material-shell-adopted-specification.md)、[採否JSON](material-shell-adoption-disposition.json)。構造全体・製品本体開始・実材料の用途別受入は未完了。

## 採用した入力と処理

等方、3D直交異方性9定数、明示全3D SPD有効則を、同じ採用済みPG-ANCFの材料点と層別U/f/K/Mへ接続する。等価入力は単層または有効層列で、C・密度・初期Greenひずみ・固定参照軸・厚さ／基準面を明示する。工学成分とMandelは変換してから評価する。法線Schurは主経路へ追加しない。

ABD単独や面密度だけでは全3D仕事／慣性が一意に定まらない反例を確認した。不足する定数や密度分布を推測して受理しない。これは任意ABD復元を必須としない既定STR-01に対応する。等価モデルの3D有効則そのものを利用者または校正処理が明示し、同じ主シェルで解析する経路を採用した。

## 証拠と監査

[事前仕様](material-shell-integration-protocol.md)に従う[追加10群](material-shell-integration-v2-result.json)が合格。等方・直交異方性・一般SPD有効則の不均等3層、非二次曲面・変厚・同じfull Vの171スカラー／513Cartesian係数を使用した。断面の二次ひずみ／速度場、初期項、基準面移動、材料軸反転、全入力経路、非線形接線、独立物理Xの慣性、界面両側回復と入力不変性を照合した。

| 指標 | 最大誤差／結果 |
|---|---:|
| 保存断面行列と独立テンソル仕事／力／質量 | {metrics['section_work_error']:.6g} |
| 同じ主シェルの接線方向微分 | {metrics['shell_tangent_error']:.6g} |
| 独立工学テンソル経路の内力 | {metrics['shell_independent_tensor_force_error']:.6g} |
| 別厚さ求積・解析物理Xによる空間慣性 | {metrics['physical_inertia_error']:.6g} |
| 材料基底反転の初期ひずみ | {metrics['initial_strain_basis_reversal_error']:.6g} |
| 入力不正の負対照 | 28件を拒否 |

[初回の統合実行](material-shell-integration-result.json)は参照評価モジュールの接続不備で停止し、6群までの結果を保持した。追加版は正しい採用参照へ接続し、同じ基準で全10群を実行した。旧結果・ソースを上書きしていない。

[既存材料点6群](material-layer-kernel-results.md)、R4の積層非線形・客観性・均質層分割・平衡復帰・質量／固有対／再開、R5の採用仕様と証拠のハッシュを確認して再利用した。主シェルの全品質試験は繰り返していない。監査は保存断面データを使った独立仕事の再計算と、C/E₀の基底変換、ソース・保存場の同一性を含む。[監査JSON](material-shell-adoption-audit.json)。

## 残る範囲

主シェルへの層別結合・等価入力・初期ひずみ・成分／軸／回復の数値仕様は1.2で閉じた。実材料の出自・有効ひずみ範囲・温度等の用途条件・実験不確かさは作業5のMaterialUseDomainと受入に対応させた。合成定数の検証を実材料の妥当性確認とはしない。

次は1.3で梁理論・DOF・軸／せん断／ねじり・大回転・断面と必要微分・質量を一組として選定し、主シェルの位置／方向・仮想仕事へ渡す境界を閉じる。一般接続・偏心慣性は1.4、荷重・減衰・回復の構造全体統合は1.5、連成時間方式は3、具体的数値依存／型／配布は4、実装開始は6で審査する。[全体ロードマップ](../roadmap.md)。
''',encoding='utf-8');sources[report.relative_to(ROOT).as_posix()]=sha(report)
    auditpath=ROOT/'docs/verification/material-shell-adoption-audit.json';links=0
    for p in [spec,report]:
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            t=(p.parent/target.split('#')[0]).resolve();require(t.exists() or t==auditpath,'missing link '+str(t));links+=1
    for p,h in sources.items():require(sha(ROOT/p)==h,'changed final source '+p)
    audit=dict(passed=True,material_1_2_complete=True,structural_review_complete=False,product_implementation_start=False,metrics=metrics,
        repeated_section_work=repeated,sources=sources,artifacts=artifacts,checked_new_links=links,seconds=time.monotonic()-tick)
    auditpath.write_text(json.dumps(audit,indent=2),encoding='utf-8');print(json.dumps(dict(metrics=metrics,sources=len(sources),artifacts=len(artifacts),links=links),indent=2))

if __name__=='__main__':main()
