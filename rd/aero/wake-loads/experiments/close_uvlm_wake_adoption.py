"""Close 2.2 while preserving all parent evidence and failed explorations."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,re
ROOT=rd_paths.ROOT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    dest=ROOT/'rd/aero/wake-loads/records/uvlm-wake-adoption-disposition.json'
    assert not dest.exists(),'immutable destination'
    audit=json.loads((ROOT/'rd/aero/kernel/records/uvlm-induction-adoption-audit.json').read_text())
    parentpath=ROOT/'rd/aero/kernel/records/uvlm-induction-adoption-closure.json';parent=json.loads(parentpath.read_text())
    assert audit['passed'] and parent['passed'] and parent['uvlm_2_1_complete']
    sources=dict(audit['sources']);artifacts=dict(audit['artifacts'])
    for group in ['immutable','status_sync_source']:sources.update(parent[group])
    sources[parentpath.relative_to(ROOT).as_posix()]=sha(parentpath)
    queue=[rd_paths.RECORDS/f'{x}-result.json' for x in ['uvlm-wake-assembly','uvlm-wake-policy']];records=[];seen=set()
    while queue:
        p=queue.pop()
        if p in seen:continue
        seen.add(p);r=json.loads(p.read_text());assert r['passed'] and not r['failures'];records+=r['records'];sources[p.relative_to(ROOT).as_posix()]=sha(p)
        for name,h in r['sources'].items():
            assert sha(ROOT/name)==h,name;sources[name]=h
            if name.endswith('-result.json'):queue.append(ROOT/name)
        if 'artifact' in r:
            ar=r['artifact'];artifacts[ar['artifact']]=ar['artifact_sha256']
    assert len(records)==16
    # Include locally imported modules recursively; installed dependencies have their own environment.
    queue=[ROOT/n for n in sources if n.endswith('.py')];seen=set()
    while queue:
        p=queue.pop()
        if p in seen:continue
        seen.add(p)
        for module in re.findall(r'^(?:from|import)\s+([A-Za-z_]\w*)',p.read_text(encoding='utf-8'),re.M):
            q=p.parent/(module+'.py')
            if q.exists():
                name=q.relative_to(ROOT).as_posix();h=sha(q)
                if name in sources:assert sources[name]==h,name
                sources[name]=h;queue.append(q)
    for stem in ['exploration','qualified','final']:
        p=ROOT/f'build/analysis/uvlm-wake-assembly-{stem}.log'
        if p.exists():artifacts[p.relative_to(ROOT).as_posix()]=sha(p)
    p=ROOT/'build/analysis/uvlm-wake-policy.log';artifacts[p.relative_to(ROOT).as_posix()]=sha(p)
    spec=ROOT/'docs/specs/aero/uvlm-wake-adopted-specification.md';assert all(f'UW-D{i:02}' in spec.read_text(encoding='utf-8') for i in range(1,9))
    decision=dict(passed=True,uvlm_2_2_complete=True,uvlm_review_complete=False,structural_review_complete=True,product_implementation_start=False,verification_groups=16,adopted_model='uvlm_wake_assembly_v2.Model',decisions=[f'UW-D{i:02}' for i in range(1,9)],next_work='2.3 independent finite-wing and unsteady references and error budgets',default_core_radius_selected=False,finite_wing_reference_qualified=False,product_time_integrator_adopted=False)
    dest.write_text(json.dumps(decision,indent=2),encoding='utf-8')
    report=ROOT/'rd/aero/wake-loads/records/uvlm-wake-adoption-results.md'
    report.write_text('''# UVLM 2.2：後流・組立・状態・全ベクトル荷重の採用

2026-11-05。**2.2を採用完了。現在地は2.3の独立精度参照・誤差予算。** [採用仕様](../uvlm-wake-adopted-specification.md)、[主12群](uvlm-wake-assembly-result.json)、[状態政策4群](uvlm-wake-policy-result.json)、[判定](uvlm-wake-adoption-disposition.json)、[監査](uvlm-wake-adoption-audit.json)。

複数の非平面揚力面、全相互誘導、陰的な後流先頭、自由／指定移流、物質ID、閉じた有限窓、全荷重、PG-ANCFへの仕事双対写像、全依存接線、試行／確定・再開の16群が合格した。採用コードは厳格な状態検査を加えたv2であり、基底の12群を変更せず追加4群で政策を照合した。

| 検証 | 結果 |
|---|---|
| 別の線・面積分による組立／循環 | 行列9.27e−16、後流基底8.11e−15、循環1.23e−15以内、元残差2.04e−16 |
| 全ベクトル荷重の独立積分 | 力・モーメント・一般化荷重の相対差2.21e−16以内 |
| 陰的循環＋全荷重の方向微分 | 相対差1.50e−12。法線・翼面速度・既知流れ・旧時刻移流を含む |
| 全自由移流の状態方向微分 | 相対差1.73e−11。自己以外の全誘導を保持 |
| 採用PGシェルの仕事・合力・モーメント | 差1.01e−15以内。PG座標の全荷重接線差7.41e−11 |
| 再開・原子性 | 両循環履歴、実刻み、モード・有限窓を保持し次状態／荷重がbitwise一致。変更済み試行、古い世代、不整合履歴等を拒否 |

相互誘導を落とすと1.32%、後流先頭を旧循環へ固定すると13.4%、自由移流の誘導を落とすと3.61%変化した。有限窓と成長後流の循環差8.04%も明示し、切捨て誤差を零とはしない。共有頂点・辺だけを許す全域面交差診断も合格した。

初回の幾何負対照はcoplanar試験点を誤って面外へ移していた。fixtureだけを訂正し、物理式や許容差は変えていない。旧ログを保存する。

Euler/BDF1は時間・履歴契約を実行検証する証拠であり、製品の時間積分方式の採用ではない。有限翼の物理荷重精度、既定半径・後流長、自由後流の時間安定性、空力全体と本体実装開始は未完了。AVLの精度参照不採用は維持する。
''',encoding='utf-8')
    for p in [spec,dest,report,Path(__file__)]:sources[p.relative_to(ROOT).as_posix()]=sha(p)
    for group in [sources,artifacts]:
        for name,h in group.items():assert sha(ROOT/name)==h,name
    ap=ROOT/'rd/aero/wake-loads/records/uvlm-wake-adoption-audit.json';ap.write_text(json.dumps(dict(passed=True,verification_groups=16,sources=sources,artifacts=artifacts),indent=2),encoding='utf-8')
    fronts=parent['mutable_status_entries']
    for name in fronts:
        p=ROOT/name;s=p.read_text(encoding='utf-8')
        s=s.replace('現在地は2.2の後流・組立・状態・荷重。','2.2の後流・組立・状態・荷重は16群の検証を経て採用完了。現在地は2.3の独立精度参照・誤差予算。')
        s=s.replace('→2.2 後流・組立・状態・荷重','→2.3 独立精度参照・誤差予算').replace('現在地2.2','現在地2.3')
        s=s.replace('**2.2 後流・組立・状態・荷重：現在地**','2.2 後流・組立・状態・荷重：採用完了')
        s=s.replace('| 2.3 独立精度参照・誤差予算 |','| **2.3 独立精度参照・誤差予算：現在地** |')
        # Preserve historical paragraphs; replace current status cells only.
        if name=='docs/specs/design-decisions.md':
            lines=s.splitlines()
            for i,line in enumerate(lines):
                if any(line.startswith('| AE-'+f'{j:02}'+' |') for j in range(1,7)):
                    cells=line.split('|');cells[-2]=' P2/P3・2.1と2.2の局所誘導・全相互組立・陰的先頭・移流・全荷重／PG写像・接線・ID／交差・状態を採用完了。16群。独立精度と既定半径／後流長は現在2.3、AE統合は2.4（[2.2仕様](uvlm-wake-adopted-specification.md)、[判定](verification/uvlm-wake-adoption-results.md)） ';lines[i]='|'.join(cells)
            s='\n'.join(lines)+'\n'
        base='docs/' if name=='README.md' else ''
        s+=f'\n2026-11-05：**2.2採用完了、現在地2.3。** [後流・全荷重仕様]({base}uvlm-wake-adopted-specification.md)、[16群の判定]({base}verification/uvlm-wake-adoption-results.md)。製品の時間方式・UVLM全体・本体実装開始は未完了。\n'
        p.write_text(s,encoding='utf-8')
    for group in [sources,artifacts]:
        for name,h in group.items():assert sha(ROOT/name)==h,name
    frozen=[spec,dest,report,ap];closure=dict(passed=True,uvlm_2_2_complete=True,uvlm_review_complete=False,product_implementation_start=False,current_work=decision['next_work'],immutable={p.relative_to(ROOT).as_posix():sha(p) for p in frozen},status_sync_source={Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))},mutable_status_entries=fronts)
    (ROOT/'rd/aero/wake-loads/records/uvlm-wake-adoption-closure.json').write_text(json.dumps(closure,indent=2),encoding='utf-8')
    print(json.dumps(dict(passed=True,sources=len(sources),artifacts=len(artifacts),verification_groups=16)))
if __name__=='__main__':main()
