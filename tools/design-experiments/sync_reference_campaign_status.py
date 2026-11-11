from pathlib import Path
import json,hashlib,re
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=ROOT/'docs/verification/uvlm-reference-campaign-audit.json';a=json.loads(ap.read_text());assert a['passed'] and not a['uvlm_2_3_complete'];dest=ROOT/'docs/verification/uvlm-reference-campaign-closure.json';assert not dest.exists()
    def verify():
        for key in ['sources','artifacts']:
            for name,h in a[key].items():assert sha(ROOT/name)==h,name
    verify();fronts=['README.md','docs/preimplementation-plan.md','docs/structural-model-adoption-review.md','docs/requirements-evidence-map.md','docs/design-decisions.md','docs/roadmap.md']
    for name in fronts:
        p=ROOT/name;s=p.read_text(encoding='utf-8');base='docs/' if name=='README.md' else ''
        if name=='docs/roadmap.md':
            lines=s.splitlines()
            for i,line in enumerate(lines):
                if line.startswith('| **2.3 独立精度参照・誤差予算：現在地** |'):lines[i]='| **2.3 独立精度参照・誤差予算：現在地** | 定常2翼の総量・20帯、k=0.2上下運動の独立参照と有限解比較、主翼＋尾翼の物理荷重比較が限定成立。最終非定常は参照不確かさ込みFz0.552%／My0.795% | [固定した結果・監査](verification/uvlm-reference-campaign-summary-2026-11-05.md)。残りは誘導抗力の共通極限、pitch／周波数等の非定常範囲、一般のコア／後流／有限解誤差規則。2.3全体は未完了 |'
            s='\n'.join(lines)+'\n'
        s+=f'\n2026-11-05・第1参照キャンペーン終了：**2.2完了、現在地2.3。** 定常2翼・20帯、非定常上下運動、主翼／尾翼の独立比較が限定成立。[成立範囲・未達・残条件・監査]({base}verification/uvlm-reference-campaign-summary-2026-11-05.md)。2.3全体・2.4・本体実装開始は未完了。\n'
        p.write_text(s,encoding='utf-8')
    p=ROOT/'docs/verification/uvlm-independent-reference-progress.md';s=p.read_text(encoding='utf-8');s=s.replace('このページは現在の整理であり、固定済みの不合格記録を変更しない。','本ページは同日の第1キャンペーン終了時点で固定する。最新の現在地はロードマップへ集約する。[固定した最終結果と監査](uvlm-reference-campaign-summary-2026-11-05.md)。')
    s=s.replace('| 帯分布の収束 | N192→384のL2差約0.215%、翼端帯の最大差約1.294%。各帯の連続極限の不確かさはまだ認定していない |','| 帯分布の収束 | N768と参照コア1e−10mまで追加し、空間・コア・比較差の最大合計指標1.7113%で2%以内。[結果](fixed-band-core-plateau-result.json)。厳密上界とはしない |')
    s=s.replace('非定常は弦M32/64と幅N192/384の方向別系列を実行し、面積・作用点の離散差と各コードの自己収束を判定する。必要な交点はM64/N192までとし、メモリ24GiB・1ケース1200秒を上限とする。[有限追加条件](unsteady-space-extension-protocol.md)。精度条件を満たした範囲だけ採用する。','非定常は方向別系列とM64/N192の交点を実施し、連続極限の独立参照を推定不確かさ約0.191%で適格化した。最後に左右対称性で縮約したM128/N192を追加し、参照不確かさ込みFz約0.552%／My約0.795%で1%以内。419秒・5.86GiBで上限内。[最終比較](symmetric-heave-final-result.json)。元の自己空間配分0.25%の合格とは区別する。')
    s=s.replace('2.3の残りは、非定常の複素荷重、固定帯の連続極限の見積り、二次の誘導抗力、複数面の独立荷重比較、既定のコア／後流長の選定規則。','2.3の残りは、pitch・換算振動数等の非定常範囲、二次の誘導抗力の共通極限、一般のコア／後流長・有限解誤差の選定規則。主翼＋尾翼の物理荷重の限定比較も成立した。')
    p.write_text(s,encoding='utf-8');verify();frozen=[ap,p,ROOT/'docs/verification/uvlm-reference-campaign-summary-2026-11-05.md'];links=0
    for name in fronts:
        q=ROOT/name
        for target in re.findall(r'\]\(([^)]+)\)',q.read_text(encoding='utf-8')):
            if 'uvlm-reference-campaign-summary' in target:assert (q.parent/target).exists();links+=1
    result=dict(passed=True,meaning='current status synchronized, evidence fixed; 2.3 remains open',uvlm_2_2_complete=True,uvlm_2_3_complete=False,uvlm_review_complete=False,product_implementation_start=False,current_work=a['current_work'],immutable={q.relative_to(ROOT).as_posix():sha(q) for q in frozen},status_sync_source={Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))},mutable_status_entries=fronts,checked_status_links=links);dest.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(dict(passed=True,status_entries=6,checked_links=links,sources=len(a['sources']),artifacts=len(a['artifacts']),uvlm_2_3_complete=False),indent=2))
if __name__=='__main__':main()
