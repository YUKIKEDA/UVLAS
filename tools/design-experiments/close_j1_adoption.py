from pathlib import Path
import json,hashlib,re
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=ROOT/'docs/verification/j1-adoption-audit.json';audit=json.loads(ap.read_text());assert audit['passed'];parent=json.loads((ROOT/'docs/verification/uvlm-integrated-adoption-closure.json').read_text());assert parent['uvlm_review_complete'];dp=ROOT/'docs/verification/j1-adoption-disposition.json';assert not dp.exists();decision=dict(passed=True,structural_review_complete=True,uvlm_review_complete=True,j1_review_complete=True,product_implementation_start=False,current_work='4 numerical foundations, concrete C++23 types and I/F, linear/eigen solvers, memory/performance and file formats',time_method='BDF2 physical velocity and circulation, implicit trapezoid newborn material points, fixed smooth segments',decisions=[f'J1-D{i:02}' for i in range(1,10)],full_practical_acceptance_complete=False);dp.write_text(json.dumps(decision,indent=2),encoding='utf-8');fronts=parent['mutable_status_entries']
    for name in fronts:
        p=ROOT/name;s=p.read_text(encoding='utf-8');base='docs/' if name=='README.md' else '';block='<!-- active-uvlm-status:start -->\n2026-11-05・最新：**構造1、UVLM2、連成・時間・安定性3の方式採用は完了。現在地は4：数値基盤・具体型・I/F・容量性能設計。** 実PGの変形平衡と全自由後流を同時に解き、全16未知量JVP、全24状態の線形化と非線形時間応答、初期境界微分、姿勢履歴、BDF2微分／全荷重の再開、枝・モード／探索制御を検証。1038ソース・642成果物を監査した。**具体数値実装、実務受入計画5、開始審査6、本体実装開始は未完了。**\n\n[J1採用仕様]('+base+'coupled-adopted-specification.md)、[終了判断・成立範囲]('+base+'verification/j1-adoption-completion.md)、[UVLM終了]('+base+'verification/uvlm-adoption-completion.md)。\n<!-- active-uvlm-status:end -->';a=s.index('<!-- active-uvlm-status:start -->');b=s.index('<!-- active-uvlm-status:end -->',a)+len('<!-- active-uvlm-status:end -->');s=s[:a]+block+s[b:]
        if name=='README.md':s=s.replace('**現在地：実装前6作業の「3. 連成・時間・安定性」。** 構造1とUVLM2は採用完了。J1共同採用・具体数値基盤・実務受入計画・実装開始審査は未完了。','**現在地：実装前6作業の「4. 数値基盤・具体型・主要I/F」。** 構造1、UVLM2、J1共同決定3は採用完了。数値基盤・実務受入計画・実装開始審査は未完了。')
        lines=s.splitlines()
        for i,line in enumerate(lines):
            if line.startswith('2026-11-05現在：構造モデル採用1とUVLM採用2'):lines[i]='2026-11-05現在：構造1、UVLM2、J1共同決定3は採用完了。現在地は4の数値基盤・具体型・主要I/F・容量性能設計。実務受入計画5、開始審査6、本体実装開始は未完了。';break
        s='\n'.join(lines)+'\n'
        if name=='docs/roadmap.md':
            a=s.index('## いま、どこにいるか');b=s.index('## 今後の計画',a);s=s[:a]+'## いま、どこにいるか\n\n**現在地は4 数値基盤・具体型・主要I/F・容量性能設計。構造1、UVLM2、J1共同決定3は方式採用完了。** [J1の終了判断](verification/j1-adoption-completion.md)、[採用仕様](coupled-adopted-specification.md)。4の具体数値解法／データ、5の実務受入計画、6の実装開始審査が残る。方式の採用は製品完成や全実務受入の完了ではない。\n\n'+s[b:];lines=s.splitlines()
            for i,line in enumerate(lines):
                if line.startswith('| **3 J1：連成・時間・安定性：現在地**'):lines[i]='| 3 J1：連成・時間・安定性 | **方式採用完了**。物理速度BDF2、新生点台形、政策別初期整合、同じ離散平衡・全履歴線形化、枝／モード・再開を[統合採用](coupled-adopted-specification.md) | 具体の疎密解法・型／形式・容量は4、実務ケースの数値／実験受入は5 | [実PG／全自由後流を含む検証と監査](verification/j1-adoption-completion.md)が成立 |'
                elif line.startswith('| 4 数値基盤・型・I/F |'):lines[i]=line.replace('| 4 数値基盤・型・I/F |','| **4 数値基盤・型・I/F：現在地** |',1)
            s='\n'.join(lines)+'\n'
        if name=='docs/design-decisions.md':
            lines=s.splitlines()
            for i,line in enumerate(lines):
                if any(line.startswith('| '+key+' |') for key in ['TR-02','TM-01','TM-02','TM-03','TM-04','TM-05','FL-01','FL-02','FL-03']):
                    cells=line.split('|');cells[-2]=' P3・J1-D01〜09の数式・状態・必要微分・拒否／精度規約を共同採用。実PG／全自由後流の変形平衡・全接線・全履歴線形化・時系列、初期微分・姿勢・再開・枝／探索を限定検証。具体数値解法／抽出器・データ実装は4、実務受入は5（[仕様](coupled-adopted-specification.md)、[終了判断](verification/j1-adoption-completion.md)） ';lines[i]='|'.join(cells)
                elif any(line.startswith('| AE-'+f'{j:02}'+' |') for j in range(1,8)):lines[i]=line.replace('AE-04/05の製品時間・初期整合はJ1共同決定3で閉じる','AE-04/05の製品時間・初期整合もJ1共同決定3で採用済み')
            s='\n'.join(lines)+'\n'
        p.write_text(s,encoding='utf-8');assert '現在地は4' in s or '現在地：実装前6作業の「4.' in s
    for name,h in audit['sources'].items():assert sha(ROOT/name)==h,name
    report=ROOT/'docs/verification/j1-adoption-completion.md'
    for target in re.findall(r'\]\(([^)]+)\)',report.read_text(encoding='utf-8')):
        if '://' not in target:assert (report.parent/target.split('#')[0]).resolve().exists(),target
    closure=dict(passed=True,structural_review_complete=True,uvlm_review_complete=True,j1_review_complete=True,product_implementation_start=False,current_work=decision['current_work'],immutable={p.relative_to(ROOT).as_posix():sha(p) for p in [ap,dp,report,ROOT/'docs/coupled-adopted-specification.md']},status_sync_source={Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))},mutable_status_entries=fronts);cp=ROOT/'docs/verification/j1-adoption-closure.json';assert not cp.exists();cp.write_text(json.dumps(closure,indent=2),encoding='utf-8');print(json.dumps(decision,indent=2))
if __name__=='__main__':main()
