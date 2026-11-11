"""Audit structural work 1 as a model-adoption gate, not implementation start."""
from pathlib import Path
import hashlib,json,re,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def require(c,m):
    if not c:raise AssertionError(m)
def main():
    start=time.monotonic();dest=ROOT/'docs/verification/structural-adoption-disposition.json';require(not dest.exists(),'immutable destination');sources={};artifacts={};parentpath=ROOT/'docs/verification/connection-adoption-closure.json';parent=json.loads(parentpath.read_text());pa=json.loads((ROOT/'docs/verification/connection-adoption-audit.json').read_text());require(parent['passed'] and parent['connection_1_4_complete'] and pa['passed'],'prerequisites')
    for group in ['sources','artifacts']:
        for name,h in pa[group].items():require(sha(ROOT/name)==h,'changed parent '+name)
        (sources if group=='sources' else artifacts).update(pa[group])
    for group in ['immutable','status_sync_source']:
        for name,h in parent[group].items():require(sha(ROOT/name)==h,'changed closure '+name);sources[name]=h
    sources[parentpath.relative_to(ROOT).as_posix()]=sha(parentpath)
    names=['structural-load-damping','structural-recovery','structural-combined','structural-history-state','structural-layer-damping'];paths=[ROOT/'docs/verification'/f'{n}-result.json' for n in names];results=[json.loads(p.read_text()) for p in paths];require(all(r['passed'] and not r['failures'] for r in results),'selected results');load,recover,combined,history,layer=results;require([len(r['records']) for r in results]==[10,4,4,5,2],'25 groups')
    for r in load['records']:
        for k,v in r.items():
            if k.endswith('_error') or k in ['rigid_motion_strain_rate']:require(v<1e-8,'load algebra '+k)
    require(load['records'][0]['symmetrization_negative_control']>.01 and load['records'][3]['closed_tangent_asymmetry']<1e-8,'open vs closed pressure');require(load['records'][7]['compressed_elastic_minimum']<0 and load['records'][7]['compressed_viscous_power']>=0,'negative stiffness separate from damping')
    cuts=recover['records'][2]['cases'];require([r['n'] for r in cuts]==[8,16,32] and max(cuts[-1]['force_error'],cuts[-1]['moment_error'])<.01,'cut accuracy');require(all(r['original_residual']<1e-8 for r in cuts),'cut original residual')
    for key in ['force_error','moment_error']:require(all(a[key]>b[key] for a,b in zip(cuts[:-1],cuts[1:])),'cut convergence')
    c=combined['records'];require(c[0]['history'][-1]['original_residual']<1e-8 and c[0]['history'][-1]['constraint']<1e-9 and c[0]['global_wrench_error']<1e-8,'actual combined equilibrium');require(c[1]['relative_error']<1e-8 and c[1]['off_equilibrium_test'],'whole residual tangent');require(max(c[2]['state_difference'],c[2]['original_residual'])<1e-8,'combined elimination');require(c[3]['reduced_mass_min']>0 and c[3]['reduced_damping_min']>=0,'matched operators')
    h=history['records'];require(h[1]['rejected']==6 and h[4]['negative_controls']==5 and h[4]['operator_difference']==0 and h[4]['bitwise_identity'],'state and histories');require(max(h[2]['response_error'],h[2]['cycle_power_error'],h[3]['tip_error'],h[3]['original_residual'])<1e-8,'independent support references')
    for r in layer['records']:
        for key,value in r.items():
            if key.endswith('_error'):require(value<1e-8,'layer damping '+key)
    queue=paths[:];seen=set()
    while queue:
        p=queue.pop()
        if p in seen:continue
        seen.add(p);result=json.loads(p.read_text());sources[p.relative_to(ROOT).as_posix()]=sha(p)
        for name,digest in result.get('sources',{}).items():
            q=ROOT/name;require(sha(q)==digest,'changed selected source '+name);sources[name]=digest
            if q.suffix=='.json':queue.append(q)
        for name,digest in result.get('artifacts',{}).items():require(sha(ROOT/name)==digest,'changed nested artifact '+name);artifacts[name]=digest
        if isinstance(result.get('artifact'),dict):
            ar=result['artifact'];require(sha(ROOT/ar['artifact'])==ar['artifact_sha256'],'snapshot');artifacts[ar['artifact']]=ar['artifact_sha256']
    a=np.load(ROOT/combined['artifact']['artifact']);require(a['M'].shape==a['K'].shape==a['C'].shape and a['J'].shape[1]==len(a['M']),'snapshot sizes');require(np.linalg.eigvalsh(a['M']).min()>0,'snapshot physical mass');require(np.linalg.norm(a['R'].T@a['R']-np.eye(3))<1e-12,'snapshot orientation')
    # Retain the rejected, inverted random coefficient perturbation and the
    # corrected admissible perturbation. Neither changes adoption tolerances.
    for name in ['structural-combined-exploration.log','structural-combined-qualified.log','structural-load-damping-exploration.log','structural-load-damping-qualified.log']:
        p=ROOT/'build/analysis'/name;require(p.exists(),'diagnostic log');artifacts[p.relative_to(ROOT).as_posix()]=sha(p)
    spec=ROOT/'docs/structural-adopted-specification.md';text=spec.read_text(encoding='utf-8');decisions=[f'SI-D{i:02d}' for i in range(1,10)];require(all(x in text for x in decisions),'decisions');require(all(f'ST-{i:02d}' in text for i in range(1,11)) and 'CV-03' in text,'ST coverage');require(all(x in text for x in ['M,qq','三階','3／4','具体数値基盤','UVLM']),'explicit handoff')
    metrics=dict(new_verification_groups=25,combined_original_residual=c[0]['history'][-1]['original_residual'],combined_constraint=c[0]['history'][-1]['constraint'],combined_wrench_error=c[0]['global_wrench_error'],whole_residual_tangent_error=c[1]['relative_error'],independent_support_series_error=h[3]['tip_error'],independent_damped_support_response_error=h[2]['response_error'],restart_operator_difference=h[4]['operator_difference'],finest_cut_force_error=cuts[-1]['force_error'],finest_cut_moment_error=cuts[-1]['moment_error'],pressure_symmetrization_defect=load['records'][0]['symmetrization_negative_control'])
    decision=dict(passed=True,structural_1_5_complete=True,structural_review_complete=True,product_implementation_start=False,main_shell_adopted=True,material_1_2_complete=True,beam_1_3_complete=True,connection_1_4_complete=True,decisions=decisions,metrics=metrics,next_work='2 UVLM model adoption',specification=spec.relative_to(ROOT).as_posix(),production_time_integration_complete=False,production_numerical_types_complete=False,real_world_acceptance_complete=False)
    dest.write_text(json.dumps(decision,indent=2),encoding='utf-8');report=ROOT/'docs/verification/structural-adoption-results.md';report.write_text(f'''# 構造モデル採用1の終了判断

2026-11-05。**1.5を完了し、1.1主シェル・1.2材料・1.3梁・1.4接続を含む構造モデル採用「1」を閉じた。** 次は2 UVLM。製品本体の実装開始条件は未達。[統合仕様v1.0](../structural-adopted-specification.md)、[判定JSON](structural-adoption-disposition.json)、[監査](structural-adoption-audit.json)、[ロードマップ](../roadmap.md)。

## 何を閉じたか

ST-01〜10/CV-03を、固定した採用式・適用／拒否・必要微分・主要処理・単位・不変条件・独立参照へ対応させた。主シェル・材料・梁・接続を探索し直さず、同じ状態に荷重・減衰・回復を統合した。選定済み構造を実装するための数値モデルを確定したという判断であり、実務構造全体の妥当性検証・製品実装完了ではない。

| 新規検証 | 結果 |
|---|---|
| [荷重・減衰10群](structural-load-damping-result.json) | 曲面の現／基準面積面力、圧力、偏心梁の世界／物体成分の力・偶力、全接線・仕事・wrench、剛体運動の内部粘性零、動く弾性／粘性支持と基準仕事 |
| [層別粘性2群](structural-layer-damping-result.json) | 異方性・異なる2層・参照上で変化する材料軸。独立な材料軸パワー、状態／速度接線、層分割の加法性 |
| [応力／切断面4群](structural-recovery-result.json) | 3D応力測度・材料軸・層片側、Nansonの面積変換、合力・基準点移動。非線形梁の切断面収支へ収束し、最細の力誤差 {metrics['finest_cut_force_error']:.4%}、偶力誤差 {metrics['finest_cut_moment_error']:.4%}、両方1%以内 |
| [実組合せ4群](structural-combined-result.json) | 採用梁―シェルに有限弾性支持・偏心質量・重力・dead荷重・追従圧力。最大最終原残差 {metrics['combined_original_residual']:.6g}、拘束 {metrics['combined_constraint']:.6g}、全物理wrench差 {metrics['combined_wrench_error']:.6g} |
| 全残差・同じM/K/C | 非平衡で全残差接線の独立微分差 {metrics['whole_residual_tangent_error']:.6g}。原KKT／独立消去が一致し、同じ平衡・拘束のM正定・C非負、追従荷重の非対称Kを保持 |
| [時系列・支持参照・状態5群](structural-history-state-result.json) | 区分線形／C1 Hermiteとイベント・片側・外挿拒否。独立な支持直列解との差 {metrics['independent_support_series_error']:.6g}、減衰支持の調和解析解との差 {metrics['independent_damped_support_response_error']:.6g}。実平衡の状態復元後、原残差／全接線の差0 |

新規25群が成立。保存された実平衡のK/M/C/J・位置／姿勢・乗数・反力も監査した。公開APIや製品数値依存は追加していない。

## 採用で落とさなかった項

開いた圧力面の接線を対称化すると、検査方向の差が {metrics['pressure_symmetrization_defect']:.2%} となった。閉じた面の体積ポテンシャルの対称性と分けて検査し、非保存成分を残した。動く支持の基準側仕事、粘性の状態接線、偏心・物体座標荷重の回転とチャート微分も保持する。圧縮平衡の負の弾性剛性を保持したまま、物理粘性の散逸非負を確認した。

最初の全残差検査では、基底係数に無尺度の乱数を加えたため試行体積が反転し、適格性検査が拒否した。許容差は変更せず、物理変形勾配を有界にした非平衡の方向で再検査した。拒否と診断ログを保存する。

## 後続へ残るもの

* 2：UVLMの幾何・渦核・自己項・荷重・後流を一組として採用し、精度参照と収束を閉じる。
* 3：構造／空力／後流の時間、双方向連成、拘束初期化、一般動的残差・高階微分、平衡周り安定性を統合する。構造単体の小系ODE照合を製品時間積分器の採用にしない。
* 4：必要精度と三階までの運動学微分、M,qq、疎構造・前処理、C++23型・所有権・寿命、外部格子・時系列／出力・チェックポイントの具体スキーマと性能。試作の密KKTやJSONを製品方式にしない。
* 5：対象翼・外板／桁／リブ・材料・接合の用途別受入、公開独立参照の条件、局所応力と面内／厚さ収束、実材料の減衰係数同定。
* 6：上記を実装パッケージへ統合し、実装開始条件を審査する。構造「1」の終了だけで本体実装を始めない。

以前の失敗や適用境界を取り消していない。学習目的や早期の部分公開を採用根拠に加えていない。
''',encoding='utf-8')
    for p in [spec,dest,report,Path(__file__)]:sources[p.relative_to(ROOT).as_posix()]=sha(p)
    # Freeze local numerical dependencies reached by imports, without hashing
    # the six deliberately mutable status entrypoints.
    pending=[ROOT/p for p in sources if p.endswith('.py')];visited=set()
    while pending:
        p=pending.pop()
        if p in visited:continue
        visited.add(p)
        for module in re.findall(r'^(?:from|import)\s+([A-Za-z_]\w*)',p.read_text(encoding='utf-8'),re.M):
            q=p.parent/(module+'.py')
            if q.exists():
                key=q.relative_to(ROOT).as_posix();digest=sha(q)
                if key in sources:require(sources[key]==digest,'changed dependency '+key)
                else:sources[key]=digest
                pending.append(q)
    auditpath=ROOT/'docs/verification/structural-adoption-audit.json';links=0
    for p in [spec,report]:
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            q=(p.parent/target.split('#')[0]).resolve();require(q.exists() or q==auditpath,'missing link '+target);links+=1
    for name,digest in sources.items():require(sha(ROOT/name)==digest,'changed final source '+name)
    auditpath.write_text(json.dumps(dict(passed=True,structural_1_5_complete=True,structural_review_complete=True,product_implementation_start=False,next_work='2 UVLM model adoption',metrics=metrics,sources=sources,artifacts=artifacts,checked_links=links,seconds=time.monotonic()-start),indent=2),encoding='utf-8');print(json.dumps(dict(passed=True,metrics=metrics,sources=len(sources),artifacts=len(artifacts)),indent=2))
if __name__=='__main__':main()
