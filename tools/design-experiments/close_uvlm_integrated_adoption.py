from pathlib import Path
import json,hashlib,re
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    auditpath=ROOT/'docs/verification/uvlm-integrated-adoption-audit.json';audit=json.loads(auditpath.read_text());assert audit['passed'];latest=json.loads((ROOT/'docs/verification/high-frequency-final-result.json').read_text());assert latest['passed'] and len(latest['records'])==12
    parent=json.loads((ROOT/'docs/verification/uvlm-reference-campaign-closure.json').read_text());fronts=parent['mutable_status_entries'];dest=ROOT/'docs/verification/uvlm-integrated-adoption-disposition.json';assert not dest.exists();decision=dict(passed=True,structural_review_complete=True,uvlm_2_1_complete=True,uvlm_2_2_complete=True,uvlm_2_3_complete=True,uvlm_2_4_complete=True,uvlm_review_complete=True,j1_review_complete=False,product_implementation_start=False,independent_harmonic_observations=12,max_combined_relative_indicator=max(r['combined_relative_indicator'] for r in latest['records']),current_work='3 joint coupling, time, initialization, constrained nonlinear mechanics, full wake state and stability',default_universal_core_selected=False,physical_product_acceptance_complete=False);dest.write_text(json.dumps(decision,indent=2),encoding='utf-8')
    current='構造モデル採用1とUVLM採用2（2.1〜2.4）は完了。現在地は3の連成・時間・安定性の共同決定。定常2翼・20帯、抗力、3周波数×上下／pitch×Fz/Myの独立参照と誤差規則を統合採用した。3の初期付加質量・線形実空力／全指定後流・非線形チャート／誕生点は限定検証済み。J1全体・具体数値基盤4・実務受入計画5・開始審査6は未完了。'
    for name in fronts:
        p=ROOT/name;s=p.read_text(encoding='utf-8');base='docs/' if name=='README.md' else '';block='<!-- active-uvlm-status:start -->\n2026-11-05・最新：**構造1とUVLM2.1〜2.4の採用は完了。現在地は3の連成・時間・安定性。** 抗力は参照不確かさ込み0.922%／0.926%、非定常12観測は全て1%以内（最大0.890%）。原本の核差と求積差を分け、初回未達は不合格のまま保持した。統合I/F、既知入力5群、大域写像4群、989ソース・625成果物の監査が成立。**J1共同採用と本体実装開始は未完了。**\n\n[UVLM終了判断]('+base+'verification/uvlm-adoption-completion.md)、[統合仕様]('+base+'uvlm-integrated-specification.md)、[連成初期化]('+base+'coupled-initialization-residual-contract.md)、[時間・安定性候補]('+base+'coupled-time-stability-specification.md)。\n<!-- active-uvlm-status:end -->'
        a=s.index('<!-- active-uvlm-status:start -->');b=s.index('<!-- active-uvlm-status:end -->',a)+len('<!-- active-uvlm-status:end -->');s=s[:a]+block+s[b:]
        # Replace the primary status paragraph, leaving dated experiment history.
        if name=='README.md':
            s=s.replace('**現在地：実装前6作業の「2. UVLM採用」→2.3 独立精度参照・誤差予算。** 構造1とUVLMの局所誘導2.1は完了。UVLM全体・連成／時間／安定性・具体数値基盤・実務受入・実装開始審査は未完了。','**現在地：実装前6作業の「3. 連成・時間・安定性」。** 構造1とUVLM2は採用完了。J1共同採用・具体数値基盤・実務受入計画・実装開始審査は未完了。')
        if name in ['README.md','docs/preimplementation-plan.md','docs/structural-model-adoption-review.md']:
            lines=s.splitlines()
            for i,line in enumerate(lines):
                if line.startswith(('2026-11-05：構造モデル採用1','現在の進行（2026-11-05）：構造モデル採用1','現在の到達点（2026-11-05）：構造モデル採用1')):lines[i]='2026-11-05現在：'+current;break
            s='\n'.join(lines)+'\n'
        if name=='docs/roadmap.md':
            a=s.index('## いま、どこにいるか');b=s.index('## 今後の計画',a);s=s[:a]+'## いま、どこにいるか\n\n**現在地は3 連成・時間・安定性。構造1とUVLM2は採用完了。** '+current+' [UVLMの終了判断](verification/uvlm-adoption-completion.md)、[構造統合仕様](structural-adopted-specification.md)。\n\n'+s[b:];s=s.replace('| **2 UVLM採用：現在地2.3** | **2.1の有限3D誘導を採用済み**。既存の限定調和応答・局所平均圧力を保持。全体採用は未完了 | **2.2後流・組立・状態・荷重→2.3独立精度参照／予算→2.4 AE統合** | AE系の採用式・履歴・荷重と独立参照、数値誤差予算・拒否条件が揃う |','| 2 UVLM採用 | **2.1〜2.4採用完了**。[統合式・主要I/F](uvlm-integrated-specification.md)、[独立参照・誤差予算・終了判断](verification/uvlm-adoption-completion.md) | 時間・高階連成は3、具体数値実装は4、実務受入は5へ引継ぎ | AE系の空間式・履歴・荷重、独立参照と誤差／拒否規則を統合済み |');s=s.replace('| 3 J1：連成・時間・安定性 |','| **3 J1：連成・時間・安定性：現在地** |')
            rows=s.splitlines()
            for i,line in enumerate(rows):
                if line.startswith('| **2.3 独立精度参照・誤差予算：現在地**'):rows[i]='| 2.3 独立精度参照・誤差予算 | **採用完了**。定常2翼・20帯、抗力2件、非定常12観測、面間相互作用、誤差規則 | [最終判定](verification/high-frequency-final-result.json)、[終了判断](verification/uvlm-adoption-completion.md)。任意ケースは別に誤差認定する |'
                elif line.startswith('| 2.4 AE統合採用 |'):rows[i]='| 2.4 AE統合採用 | **採用完了**。AE-01〜07の式・所有権・必要微分・主要I/F、入力5群・大域写像4群 | [統合仕様](uvlm-integrated-specification.md)。3の共同決定へ引継ぎ |'
            s='\n'.join(rows)+'\n'
        if name=='docs/design-decisions.md':
            lines=s.splitlines()
            for i,line in enumerate(lines):
                if any(line.startswith('| AE-'+f'{j:02}'+' |') for j in range(1,8)):
                    cells=line.split('|');cells[-2]=' P3・空間方式／入力／主要契約を2.1〜2.4で統合採用。定常・分布・抗力・非定常12観測と誤差規則、全相互作用・全荷重／PG写像・状態／再開を監査。AE-04/05の製品時間・初期整合はJ1共同決定3で閉じる（[統合仕様](uvlm-integrated-specification.md)、[終了判断](verification/uvlm-adoption-completion.md)） ';lines[i]='|'.join(cells)
            s='\n'.join(lines)+'\n'
        p.write_text(s,encoding='utf-8');assert '現在地は3' in s or '現在地：実装前6作業の「3.' in s
    # Read immutable sources once more after the status-only changes.
    for name,h in audit['sources'].items():assert sha(ROOT/name)==h,name
    report=ROOT/'docs/verification/uvlm-adoption-completion.md'
    for target in re.findall(r'\]\(([^)]+)\)',report.read_text(encoding='utf-8')):
        if '://' not in target:assert (report.parent/target.split('#')[0]).resolve().exists(),target
    closure=dict(passed=True,uvlm_review_complete=True,uvlm_2_3_complete=True,uvlm_2_4_complete=True,j1_review_complete=False,product_implementation_start=False,current_work=decision['current_work'],immutable={p.relative_to(ROOT).as_posix():sha(p) for p in [auditpath,dest,report,ROOT/'docs/uvlm-integrated-specification.md']},status_sync_source={Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))},mutable_status_entries=fronts);cp=ROOT/'docs/verification/uvlm-integrated-adoption-closure.json';assert not cp.exists();cp.write_text(json.dumps(closure,indent=2),encoding='utf-8');print(json.dumps(decision,indent=2))
if __name__=='__main__':main()
