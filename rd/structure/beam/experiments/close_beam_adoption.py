"""Immutable beam adoption audit, including earlier failures and precision gates."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,re,time
import numpy as np
import beam_se3_v2 as b
from check_beam_se3_physical_v4 import section
ROOT=rd_paths.ROOT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(c,m):
    if not c:raise AssertionError(m)
def main():
    start=time.monotonic();dest=ROOT/'rd/structure/beam/records/beam-adoption-disposition.json';require(not dest.exists(),'immutable destination')
    previous=ROOT/'rd/structure/material/records/material-shell-adoption-audit.json';pa=json.loads(previous.read_text());pc=ROOT/'rd/structure/material/records/material-shell-adoption-closure.json';parent=json.loads(pc.read_text());require(parent['material_1_2_complete'] and pa['passed'],'material prerequisites')
    sources={};artifacts={}
    for name,h in parent['immutable'].items():require(sha(ROOT/name)==h,'changed material closure '+name);sources[name]=h
    for label in ['sources','artifacts']:
        for name,h in pa[label].items():require(sha(ROOT/name)==h,'changed parent evidence '+name)
        (sources if label=='sources' else artifacts).update(pa[label])
    sources[pc.relative_to(ROOT).as_posix()]=sha(pc)
    names=['beam-se3-core-v2','beam-se3-physical-v4','beam-section-interface','beam-guards-precision'];reports={name:json.loads((rd_paths.RECORDS/f'{name}-result.json').read_text()) for name in names}
    require(all(r['passed'] and not r.get('failures') for r in reports.values()),'selected beam results')
    core=reports[names[0]];physical=reports[names[1]];interface=reports[names[2]];guards=reports[names[3]]
    require([len(r['records']) for r in reports.values()]==[9,6,3,2],'beam coverage')
    c={r['case']:r for r in core['records']};p={r['case']:r for r in physical['records']}
    for r in core['records']+interface['records']:
        if 'metrics' in r:require(max(r['metrics'].values())<=1e-8,'directional algebra')
    require(c['free_element_rigid_six']['root_ranks']==[6]*3 and c['free_element_rigid_six']['nullity']==6,'six rigid modes')
    require(c['consistent_eccentric_mass_and_physical_inertia']['mass_min']>0 and c['consistent_eccentric_mass_and_physical_inertia']['spatial_inertia_error']<=1e-9,'physical mass')
    require(c['configuration_mass_and_Christoffel']['mass_derivative_error']<=1e-7 and c['configuration_mass_and_Christoffel']['power_identity_error']<=1e-10,'gyro')
    bends=p['large_bending_two_axes_three_slenderness']['cases'];require(len(bends)==12 and all(max(s['strain_error'],s['wrench_dual_error'])<=1e-8 and s['circle_error']<=1e-10 for s in bends),'pure bending')
    require(sum(s['precision_digits']==80 for s in bends)==4,'high precision poses')
    linear=p['Timoshenko_tip_force_convergence']['cases'];require(len(linear)==12 and all(s['residual']<=1e-8 for s in linear),'linear residual')
    for h in [.1,.01,.0001]:
        seq=[s['error'] for s in linear if s['h']==h];require(seq[-1]<=.01 and all(x>=y for x,y in zip(seq[:-1],seq[1:])),'linear convergence')
    for n in [4,8,16,32]:
        errors=[s['error'] for s in linear if s['n']==n];require(max(errors)/min(errors)<=2,'slenderness amplification')
    nonlinear=p['nonlinear_tip_force_vs_continuum_BVP']['cases'];require(len(nonlinear)==4 and all(s['residual']<=1e-8 for s in nonlinear),'nonlinear force residual')
    require(max(nonlinear[-1][k] for k in ['position_error','angle_error'])<=.01,'nonlinear independent BVP accuracy')
    modes=p['consistent_modes_vs_analytic_and_continuum']['cases'];require(len(modes)==9,'modal coverage')
    for s in modes:
        require(max(s['backward'],s['mass_orthogonality'])<=1e-8,'modal algebra')
        if s['n']==64:require(max(s['frequency_errors'])<=.004,'modal accuracy')
    torsion=p['Saint_Venant_torsion_vs_2D_Poisson'];require(torsion['cases'][-1]['relative_error']<=.01 and torsion['polar_substitution_defect']>.1,'independent cross-section')
    require(interface['records'][2]['force_moment_wrench_error']<=1e-10 and interface['records'][2]['virtual_work_error']<=1e-10,'shell point/frame transfer')
    require(len(guards['records'][0]['negative_controls'])==7,'geometry negative controls')
    require(len(guards['records'][1]['cases'])==3 and all(s['retained_residual']<=1e-8 and not s['projection_qualified'] for s in guards['records'][1]['cases']),'precision retention and projection rejection')
    # Traverse source-bearing result JSONs, preserving failed variants as parents.
    queue=[rd_paths.RECORDS/f'{name}-result.json' for name in names];seen=set()
    while queue:
        path=queue.pop()
        if path in seen:continue
        seen.add(path);r=json.loads(path.read_text());sources[path.relative_to(ROOT).as_posix()]=sha(path)
        for name,h in r.get('sources',{}).items():
            q=ROOT/name;require(sha(q)==h,'changed beam source '+name);sources[name]=h
            if q.suffix=='.json':queue.append(q)
        if isinstance(r.get('artifact'),dict):
            ar=r['artifact'];require(sha(ROOT/ar['artifact'])==ar['artifact_sha256'],'beam snapshot');artifacts[ar['artifact']]=ar['artifact_sha256']
    snapshot=core['artifact'];a=np.load(ROOT/snapshot['artifact']);require(a['tangent'].shape==(12,12) and a['mass_derivative'].shape==(12,12,12),'saved operator dimensions')
    require(np.linalg.norm(a['tangent']-a['tangent'].T)<=1e-10*np.linalg.norm(a['tangent']) and np.linalg.eigvalsh(a['mass']).min()>0,'saved tangent and mass')
    # Material-coordinate reparameterization: all data, not chord length alone.
    alpha=1.7;ref=a['reference'];e0=a['e0'];K=a['K_section'];L=float(a['L']);r=a['r'];R=a['R'];x=b.evaluate(r,R,L,K,ref,e0);y=b.evaluate(r,R,L/alpha,K/alpha,alpha*ref,alpha*e0)
    coordinate_error=max(abs(x['U']-y['U'])/abs(x['U']),np.linalg.norm(x['f']-y['f'])/np.linalg.norm(x['f']),np.linalg.norm(x['K']-y['K'])/np.linalg.norm(x['K']));require(coordinate_error<=1e-10,'reference metric conversion')
    # Supported compressive equilibrium retains actual negative stiffness.
    Ks,_,_=section(.1);straight=np.array([[0.,0,0],[.99,0,0]]);eye=np.array([np.eye(3)]*2);state=b.evaluate(straight,eye,1.,Ks,np.array([1.,0,0,0,0,0]));columns=[7,8,9,10,11];values=np.linalg.eigvalsh(state['K'][np.ix_(columns,columns)]);require(values.min()<0,'negative physical stiffness retained')
    supplemental=dict(reference_coordinate_energy_force_tangent_error=float(coordinate_error),compression_negative_eigenvalues=values[values<0].tolist())
    spec=ROOT/'docs/specs/structure/beam-adopted-specification.md';text=spec.read_text();decisions=[f'BM-D{i:02d}' for i in range(1,10)];require(all(s in text for s in decisions),'adoption decisions')
    for path in [spec,Path(__file__),ROOT/'docs/requirements/product-requirements.md',ROOT/'docs/specs/structure/main-shell-adopted-specification.md',ROOT/'docs/specs/structure/material-shell-adopted-specification.md']:
        sources[path.relative_to(ROOT).as_posix()]=sha(path)
    metrics=dict(verification_groups=20,input_negative_controls=25,max_pure_bending_wrench_dual_error=max(s['wrench_dual_error'] for s in bends),
        finest_linear_tip_error=max(s['error'] for s in linear if s['n']==32),finest_nonlinear_tip_error=nonlinear[-1]['position_error'],
        finest_modal_error=max(max(s['frequency_errors']) for s in modes if s['n']==64),torsion_section_error=torsion['cases'][-1]['relative_error'],
        physical_inertia_error=c['consistent_eccentric_mass_and_physical_inertia']['spatial_inertia_error'],max_rejected_binary64_projection=max(s['binary64_projection_residual'] for s in guards['records'][1]['cases']))
    decision=dict(beam_1_3_complete=True,adopted_model='two-node SE(3) geometrically exact shear-deformable beam',decisions=decisions,
        main_shell_adopted=True,material_1_2_complete=True,structural_review_complete=False,product_implementation_start=False,metrics=metrics,
        constrained_warping_supported=False,local_stress_recovery='only with qualified sectional recovery data',
        precision_requirement='retain sufficient pose/state precision and qualify output projection; concrete production types in work 4',
        next_work='1.4 connections supports mass and eccentric inertia',specification=spec.relative_to(ROOT).as_posix(),real_material_validation=False)
    dest.write_text(json.dumps(decision,indent=2),encoding='utf-8');sources[dest.relative_to(ROOT).as_posix()]=sha(dest)
    report=ROOT/'rd/structure/beam/records/beam-adoption-results.md'
    report.write_text(f'''# 梁1.3完了：SE(3)梁の採用

2026-11-05。**SE(3)幾何学的厳密せん断変形梁を採用し、1.3の理論・自由度・断面・質量・シェル受渡しを固定した。** [採用仕様v1.0](../beam-adopted-specification.md)、[採否JSON](beam-adoption-disposition.json)。次は1.4接続・支持・質量・偏心慣性。構造全体・本体実装開始は未完了。

## 比較と成立根拠

[候補比較・事前条件](beam-se3-adoption-protocol.md)で、せん断を含まない共回転Euler–Bernoulliと、断面3D補間・対策を同時に選ぶ必要があるANCF梁を主候補とせず、位置と回転を連成補間するSE(3)梁を選んだ。梁までANCF座標に限定する要求ではない。公刊補間の出典は [Sonneville–Cardona–Brüls 2014](https://orbi.uliege.be/handle/2268/159471)。外部コードはコピーせず、今回の増分・変分・微分・慣性を独立に実装した。

[基本9群](beam-se3-core-v2-result.json)、[物理問題6群](beam-se3-physical-v4-result.json)、[断面／受渡し3群](beam-section-interface-result.json)、[正則性／精度2群](beam-guards-precision-result.json)の計20群が合格。入力負対照は25件。さらに原点・参照物質座標の変換と圧縮平衡の負固有値保持を監査した。

| 審査対象 | 結果 |
|---|---|
| 大回転・薄肉純曲げ | 2軸、厚さ比.1/.01/.0001、1要素1.8radと分割した全体4.5radの12配置。最大断面双対誤差 {metrics['max_pure_bending_wrench_dual_error']:.6g} |
| 線形先端力 | n4/8/16/32、3厚さの12問題。最細変位誤差最大 {metrics['finest_linear_tip_error']:.6g}、細分化減少と薄さ増幅2倍以内 |
| 大変形先端力 | 原Newton残差1e−8以内。独立な連続Cosserat境界値ODEへ収束し、最細先端位置誤差 {metrics['finest_nonlinear_tip_error']:.6g} |
| 整合M/Kの固有対 | 軸・ねじりの固定―自由解析解と曲げの独立Timoshenko伝達ODEへ比較。最細対象モード誤差最大 {metrics['finest_modal_error']:.6g}、0.4%以内 |
| 断面ねじり | 長方形の独立2D Poisson問題へ収束、最細差 {metrics['torsion_section_error']:.6g}。極二次モーメントによる代用を負対照で検出 |
| 慣性・受渡し | 偏心断面の物理空間慣性、M(q)とChristoffel、点／フレームのJ/H・転置仕事・偏心偶力が成立 |

いずれも合成弾性断面と規定した連続体モデルの検証である。大きな回転の検査を、その曲率・断面寸法での実材料の有限ひずみ妥当性とはしない。

## 不合格を変えた数値条件

[初期float64版](beam-se3-physical-result.json)と[等価式版](beam-se3-physical-v2-result.json)では、極細梁の純曲げで小さい不要ひずみがモーメントへ増幅された。[次版](beam-se3-physical-v3-result.json)は純曲げの状態・評価精度を保持したが、線形求解の原残差で停止した。全て保存し、上書きしていない。

採用版は同じSE(3)補間・断面則・閾値を使い、必要な姿勢は生成時から精度を保ち、U/f/Kを高精度で評価する。線形系も同じ剛性平方根で解き、原残差を評価するまで状態を保持する。[精度条件](beam-se3-precision-addendum.md)。高精度状態をfloat64へ投影した3ケースでは原残差が1e−8を超え、最大 {metrics['max_rejected_binary64_projection']:.6g}。その投影を次の有効状態へ使わない規約を固定した。閾値を緩めて採用したものではない。

Decimal80／Arb256は設計試作の証拠で、製品依存や全処理の完成結果ではない。具体的な精度型・状態保存・混合精度・投影後保証・疎求解と性能は4、時間残差と必要な高階微分は3／4の開始条件へ渡す。必要精度を満たせない実装はその条件で停止する。

## 適用境界と次作業

主経路は自由Saint-Venant反りの断面効果をK₆に含む6自由度断面であり、拘束反り・バイモーメント・断面変形が必要な入力は拒否する。局所応力は適格な断面回復データがある場合だけ返す。これらを梁の対象範囲へ黙って含めない。参照／現在の断面管、相対Log、物理的質量、成分と単位を検査する。

シェルへの物理点・フレーム・一次／二次微分・転置力を固定したが、実際の梁―シェル・折れ・交差・剛体接合の荷重経路・反力・固有特性は次の1.4で比較する。入力STR-03、偏心慣性STR-04、材料・出力要求は維持する。[全体ロードマップ](../roadmap.md)、[監査](beam-adoption-audit.json)。
''',encoding='utf-8');sources[report.relative_to(ROOT).as_posix()]=sha(report)
    auditpath=ROOT/'rd/structure/beam/records/beam-adoption-audit.json';links=0
    for path in [spec,report]:
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text()):
            if '://' in target or target.startswith('#'):continue
            q=(path.parent/target.split('#')[0]).resolve();require(q.exists() or q==auditpath,'missing link '+str(q));links+=1
    for p,h in sources.items():require(sha(ROOT/p)==h,'changed final evidence '+p)
    audit=dict(passed=True,beam_1_3_complete=True,structural_review_complete=False,product_implementation_start=False,metrics=metrics,supplemental=supplemental,
        sources=sources,artifacts=artifacts,checked_new_links=links,external_sources=[dict(url='https://orbi.uliege.be/handle/2268/159471',role='SE(3) interpolation theory',year=2014),dict(url='https://api.projectchrono.org/9.0.0/manual_fea_elements.html',role='comparison ANCF and Euler'),dict(url='https://ansyshelp.ansys.com/public/views/secured/corp/v261/en/ans_elem/Hlp_E_BEAM188.html',role='restrained warping distinction')],seconds=time.monotonic()-start)
    auditpath.write_text(json.dumps(audit,indent=2),encoding='utf-8');print(json.dumps(dict(metrics=metrics,supplemental=supplemental,sources=len(sources),artifacts=len(artifacts)),indent=2))
if __name__=='__main__':main()
