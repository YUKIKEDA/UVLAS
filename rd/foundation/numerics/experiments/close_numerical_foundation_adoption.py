import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json,hashlib,re
ROOT=rd_paths.ROOT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    audit_path=ROOT/'rd/foundation/numerics/records/numerical-foundation-adoption-audit.json';audit=json.loads(audit_path.read_text());assert audit['passed'];parent=json.loads((ROOT/'rd/coupled/j1/records/j1-adoption-closure.json').read_text());fronts=parent['mutable_status_entries']
    for name in fronts:
        p=ROOT/name;s=p.read_text(encoding='utf-8');base='docs/' if name=='README.md' else ''
        block='<!-- active-uvlm-status:start -->\n2026-11-06・最新：**構造1、UVLM2、連成・時間・安定性3、数値基盤・具体型・主要I/F4は採用設計完了。現在地は5：実務ケース・公開参照・受入／性能計画。** 4では両OSの数値4群＋Python接続1群ずつが成立し、1058ソース・671成果物を監査。5,000密／100k合成疎の容量を確認したが、実PG全系・自由後流・探索全体の製品性能は未測定。**5の終了、6の開始審査、本体実装開始、実務受入は未完了。**\n\n[4の採用仕様]('+base+'numerical-foundation-adoption-specification.md)、[入出力仕様]('+base+'io-execution-adoption-specification.md)、[終了判断・成立範囲]('+base+'verification/numerical-foundation-adoption-completion.md)。\n<!-- active-uvlm-status:end -->'
        s=re.sub(r'<!-- active-uvlm-status:start -->[\s\S]*?<!-- active-uvlm-status:end -->',lambda _:block,s,count=1)
        if name=='README.md':s=s.replace('**現在地：実装前6作業の「4. 数値基盤・具体型・主要I/F」。** 構造1、UVLM2、J1共同決定3は採用完了。数値基盤・実務受入計画・実装開始審査は未完了。','**現在地：実装前6作業の「5. 実務ケース・受入計画」。** 1〜4の方式・具体契約は採用完了。実務受入計画・全体開始審査・本体実装開始は未完了。')
        if name=='rd/planning/preimplementation-history.md':
            a=s.index('## いま、どこにいるか');b=s.index('## 今後の計画',a)
            s=s[:a]+'## いま、どこにいるか\n\n**現在地は5 実務ケース・公開参照・受入／性能計画。** 1構造、2UVLM、3J1、4数値基盤・具体型・I/Fの採用設計を終了した。[4の終了判断](verification/numerical-foundation-adoption-completion.md)。残るのは、5で実装後に何をどの入力・参照・許容差で受け入れるかを固定し、6で全体の実装開始条件を審査すること。製品性能・実務受入の実測完了は本体実装後であり、開始条件と混同しない。\n\n'+s[b:]
            s=s.replace('C[2 UVLM採用：現在地]','C[2 UVLM採用：完了]').replace('D[3 J1 連成・時間・安定性]','D[3 J1採用：完了]').replace('E[4 数値基盤・型・I/F]','E[4 数値基盤・型・I/F：完了]').replace('F[5 実務ケース・受入根拠]','F[5 実務ケース・受入計画：現在地]')
            s=s.replace('構造1.1〜1.5の採用は終了。現在は2 UVLMを進め、その後3連成・時間・安定性、4数値基盤・主要I/F、5実務ケース・受入計画、6実装開始審査へ進む。','1構造、2UVLM、3連成・時間・安定性、4数値基盤・主要I/Fの採用設計は終了。現在は5実務ケース・受入計画、その後6実装開始審査へ進む。')
            s=s.replace('構造モデル採用1を終了し、次は2 UVLMの方式・精度参照を閉じる。','この時点で構造モデル採用1を終了した。後続2〜4も現在は採用完了。')
            s=s.replace('現在の2 UVLMで構造の候補探索を再開しない。','採用済み1〜4の候補探索を根拠なく再開しない。')
            lines=s.splitlines()
            for i,line in enumerate(lines):
                if line.startswith('| **4 数値基盤・型・I/F：現在地** |'):lines[i]='| 4 数値基盤・型・I/F | **採用設計完了**。N4-01〜10／IO4-01〜07、両OSの数値・binding試作、容量を監査 | 製品への移植と全消費者・実graph・配布検証は実装時のgate | [終了判断](verification/numerical-foundation-adoption-completion.md) |'
                elif line.startswith('| 5 実務'):lines[i]=line.replace('| 5 実務','| **5 実務',1).replace(' |','：現在地** |',1)
            s='\n'.join(lines)+'\n'
        p.write_text(s,encoding='utf-8')
    for name,h in audit['sources'].items():assert sha(ROOT/name)==h,name
    decision=dict(passed=True,structural_review_complete=True,uvlm_review_complete=True,j1_review_complete=True,numerical_foundation_review_complete=True,product_implementation_start=False,practical_acceptance_complete=False,current_work=audit['current_work'])
    dest=ROOT/'rd/foundation/numerics/records/numerical-foundation-adoption-disposition.json';assert not dest.exists();dest.write_text(json.dumps(decision,indent=2),encoding='utf-8')
    closure=dict(**decision,immutable={p.relative_to(ROOT).as_posix():sha(p) for p in [audit_path,dest,ROOT/'rd/foundation/numerics/records/numerical-foundation-adoption-completion.md']},status_sync_source={Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))},mutable_status_entries=fronts)
    dest=ROOT/'rd/foundation/numerics/records/numerical-foundation-adoption-closure.json';assert not dest.exists();dest.write_text(json.dumps(closure,indent=2),encoding='utf-8');print(json.dumps(decision))
if __name__=='__main__':main()
