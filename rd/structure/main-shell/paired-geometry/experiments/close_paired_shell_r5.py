"""R5 adoption: verify frozen R4 + swept taper correction, publish decision."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,re,time
import numpy as np
ROOT=rd_paths.ROOT
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def check(condition,message):
    if not condition:raise AssertionError(message)
def pairs(x):
    if isinstance(x,dict):
        for name in ['artifact','parent_artifact']:
            if name in x and name+'_sha256' in x:yield x[name],x[name+'_sha256']
        for v in x.values():yield from pairs(v)
    elif isinstance(x,list):
        for v in x:yield from pairs(v)

def main():
    tick=time.monotonic();dest=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-shell-r5-disposition.json';check(not dest.exists(),'immutable destination')
    r4path=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-reference-r4-audit.json';r4=json.loads(r4path.read_text())
    r4close=json.loads((ROOT/'rd/structure/main-shell/paired-geometry/records/paired-reference-r4-closure.json').read_text());check(r4['passed'] and r4close['R4_complete'],'R4 prerequisites')
    manifest={};artifacts={}
    for k,v in r4close['immutable'].items():check(sha(ROOT/k)==v,'changed R4 closure '+k);manifest[k]=v
    for label in ['sources','artifacts']:
        for k,v in r4[label].items():check(sha(ROOT/k)==v,'changed R4 '+k)
        (manifest if label=='sources' else artifacts).update(r4[label])
    previous=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-shell-r5-domain-result.json';parent=json.loads(previous.read_text())
    path=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-shell-r5-domain-v2-result.json';r=json.loads(path.read_text());check(r['passed'] and len(r['records'])==4,'domain result')
    check(not parent['passed'],'retain initial load integration failure')
    for report_path,report in [(previous,parent),(path,r)]:
        manifest[report_path.relative_to(ROOT).as_posix()]=sha(report_path)
        for p,h in report['sources'].items():check(sha(ROOT/p)==h,'changed source '+p);manifest[p]=h
        for p,h in pairs(report):check(sha(ROOT/p)==h,'changed artifact '+p);artifacts[p]=h
    for record in r['records']:
        ns=record['ns'];check(record['stiffness_ranks']==[3*ns]*3 and record['mass_root_ranks']==[ns]*3,'domain rank')
        check(record['rcond']>=1e-10 and record['mass_min']>1e-10 and record['cold_quality_seconds']<=300,'domain numerical conditions')
        check(max(record['observation']['stationarity'])<=1e-8 and max(record['physical']['stationarity'])<=1e-8,'original stationarity')
        check(max(record['physical']['error'])<=.01 and max(record['physical']['position'])<=.01,'physical accuracy')
        check(record['physical_vs_weak']<=1e-8 and record['work_error']<=.01 and np.linalg.norm(record['loads']['resultant'])<=1e-10,'physical load gate')
        check(np.max(abs(np.array(record['footprint_moments'])-[.8,11/30,131/300]))<=1e-13,'exact footprint moments')
        # Recompute solve residual and independent-load solution difference from saved matrices.
        a=np.load(ROOT/record['parent_artifact']);b=np.load(ROOT/record['artifact']);T=a['root'];scale=a['scale'];coef=b['physical_coef'];weakcoef=b['weak_coef']
        rhs=scale[:,None]*b['force'].reshape(-1,1)/record['h']**3
        backward=float(np.linalg.norm(T.T@(T@(coef/scale[:,None]))-rhs)/np.linalg.norm(rhs))
        error=float(np.linalg.norm(T@((coef-weakcoef)/scale[:,None]))/np.linalg.norm(T@(weakcoef/scale[:,None])))
        check(max(backward,error)<=1e-8,'saved physical force solve')
        check(abs(error-record['physical_vs_weak'])<=1e-14,'saved load metric reproduction')
    for n in [1,2]:
        a,b=[s for s in r['records'] if s['n']==n]
        for key in ['error','position']:check(np.all(np.maximum(b['observation'][key],1e-12)/np.maximum(a['observation'][key],1e-12)<=2),'thin amplification')
    for h in [.01,.0001]:
        a,b=[s for s in r['records'] if s['h']==h]
        for key in ['error','position']:
            values=np.array(b['observation'][key]);check(np.all(values<=.01) and np.all(values<=np.maximum(a['observation'][key],1e-12)),'finest accuracy and decrease')
    spec=ROOT/'docs/specs/structure/main-shell-adopted-specification.md';handoff=ROOT/'rd/planning/structural-integration-after-shell.md'
    text=spec.read_text(encoding='utf-8');decisions=[f'SH-D{i:02d}' for i in range(1,9)]
    check(all(d in text for d in decisions),'technical decision coverage')
    check(all(s in text for s in ['1/L_ref','T_aL_a','3i+','C_M','KKT','拒否','参照識別子']),'specification coverage')
    for p in [spec,handoff,Path(__file__),ROOT/'docs/requirements/product-requirements.md',ROOT/'docs/requirements/acceptance-criteria.md',ROOT/'docs/specs/structure/material-layer-adoption-contract.md']:
        manifest[p.relative_to(ROOT).as_posix()]=sha(p)
    metrics=dict(additional_quality_configurations=4,additional_quality_rhs=12,additional_physical_problems=4,
        finest_quality_energy=max(max(s['observation']['error']) for s in r['records'] if s['n']==2),
        finest_quality_position=max(max(s['observation']['position']) for s in r['records'] if s['n']==2),
        max_thin_amplification=max(max(c['values']) for c in r['comparisons'] if c['kind']=='thin_amplification'),
        max_physical_vs_weak=max(s['physical_vs_weak'] for s in r['records']),max_physical_stationarity=max(max(s['physical']['stationarity']) for s in r['records']))
    disposition=dict(R5_complete=True,main_shell_adopted=True,adopted_model='reference-coefficient PG-ANCF',specification='docs/specs/structure/main-shell-adopted-specification.md',
        decisions=decisions,scope='main-shell numerical formulation and internal contracts; product acceptance remains required',
        structural_review_complete=False,product_implementation_start=False,next_work='1.2 material/layer/equivalent-input integration',metrics=metrics,
        evidence=['rd/structure/main-shell/paired-geometry/records/paired-reference-r4-audit.json','rd/structure/main-shell/paired-geometry/records/paired-reference-r4-closure.json','rd/structure/main-shell/paired-geometry/records/paired-shell-r5-domain-v2-result.json'],
        required_downstream={'1.2':'material and equivalent input','1.3':'beam adoption','1.4':'connections supports and eccentric inertia','1.5':'loads damping recovery and structural integration',
            '2':'UVLM adoption','3':'J1 time and stability','4':'sparse numerical dependencies precision and concrete types','5':'practical acceptance cases geometry error material validity and connected performance','6':'implementation start review'},
        unchanged_product_requirements=sha(ROOT/'docs/requirements/product-requirements.md'))
    dest.write_text(json.dumps(disposition,indent=2),encoding='utf-8');manifest[dest.relative_to(ROOT).as_posix()]=sha(dest)
    report=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-shell-r5-results.md'
    report.write_text(f'''# R5完了：主シェルの採用と仕様固定

2026-11-05。**参照係数PG-ANCFを主シェルとして採用し、R5と構造1.1の採用判断を完了する。** [採用仕様v1.0](../main-shell-adopted-specification.md)、[採否JSON](paired-shell-r5-disposition.json)。次は[1.2材料・積層・等価入力の統合](../structural-integration-after-shell.md)。構造1全体・製品本体の実装開始は未完了。

## 判断の根拠

[R4の9群・28配置84右辺](paired-reference-r4-results.md)を再実行せず、全ソース・保存場の同一性を確認した。R5では製品要求GEO-01との対応のため、後退・テーパの台形領域を追加した。同じ曲面・変厚・位置空間でn1/2・2厚さ、独立曲率a3/5/7の4配置12右辺、および独立な体積力／全6面力の4物理問題を検査した。[追加仕様](paired-shell-r5-domain-protocol.md)、[採用する追加結果](paired-shell-r5-domain-v2-result.json)。

追加の最細品質誤差はエネルギー {metrics['finest_quality_energy']:.6g}、位置 {metrics['finest_quality_position']:.6g}。薄肉増幅最大 {metrics['max_thin_amplification']:.6g} 倍。物理荷重解と独立弱積分解との差は最大 {metrics['max_physical_vs_weak']:.6g}、原停留は最大 {metrics['max_physical_stationarity']:.6g}。1%・2倍・原残差1e−8等の既定基準を通過した。

初回の[台形荷重試験](paired-shell-r5-domain-result.json)は、幅が変わった先端側面のGauss区間を背景格子の実際のt境界で区切っておらず、物理荷重の原停留が未達だった。交点分割を修正し、保存した同じN・剛性因子・尺度をそのまま使って独立荷重を再積分・再求解した。旧結果とソースは保持し、空間・材料・許容差を変えて合格にしたものではない。この条件をSH-D07として採用仕様へ反映した。

## 今回固定したもの

SH-D01〜08として、参照直接持上げ、r/d/cとP7、一定Cartesian位置空間・整合質量、有理係数と参照条件付け、全3D積層弾性、全位置トレースと転置力、有界パッチと疎接続、非矩形領域の積分、同じ平衡のM/Kを固定した。採用仕様には単位、並び、全必要微分、内力／材料＋幾何接線、求積、参照／基底識別子、所有者とworker、失敗保証、回復、状態、費用を対応させた。

採用するのは規定した数値モデルと内部契約である。旧膜Ritz固定案の不採用を取り消さない。公開API・数値依存・CLI配布・製品時間方式を今回固定したとはしない。任意のトリム配置、全実務材料、交差・梁接続、外部格子変換、全結合10万自由度性能、連成・時間・フラッタの受入は[担当と終了条件](../structural-integration-after-shell.md)へ明示した。STR-01等価特性、STR-02外板／桁／リブ、STR-03梁、STR-04偏心慣性等の必須要求は維持する。

## 更新と追跡

README、ロードマップ、実装前計画、構造審査、要求証拠表、ST-01/04の現在地を同期する。R4の固定記録は変更しない。新しい反例は採用決定IDと対応させ、根拠なく候補探索や完了したR4を再開しない。[R5監査](paired-shell-r5-audit.json)、[全体ロードマップ](../roadmap.md)。
''',encoding='utf-8');manifest[report.relative_to(ROOT).as_posix()]=sha(report)
    # Check new links before publication. The audit below is the only pending target.
    audit_path=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-shell-r5-audit.json';links=0
    for p in [spec,handoff,report]:
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            link=(p.parent/target.split('#')[0]).resolve();check(link.exists() or link==audit_path,'broken link '+str(link));links+=1
    for p,h in manifest.items():check(sha(ROOT/p)==h,'changed final evidence '+p)
    audit=dict(passed=True,R5_complete=True,main_shell_adopted=True,structural_review_complete=False,product_implementation_start=False,
        metrics=metrics,sources=manifest,artifacts=artifacts,checked_new_links=links,seconds=time.monotonic()-tick)
    audit_path.write_text(json.dumps(audit,indent=2),encoding='utf-8');print(json.dumps(dict(metrics=metrics,sources=len(manifest),artifacts=len(artifacts),links=links),indent=2))

if __name__=='__main__':main()
