"""Update the six mutable current-status entry points after immutable R5 audit."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,re
ROOT=rd_paths.ROOT
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    audit_path=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-shell-r5-audit.json';a=json.loads(audit_path.read_text());assert a['passed'] and a['main_shell_adopted']
    for group in ['sources','artifacts']:
        for name,h in a[group].items():assert sha(ROOT/name)==h,name
    common='R5を完了し、参照係数PG-ANCFを主シェルとして採用。適用／拒否条件、位置空間・必要微分、全3D材料／質量、全トレース、求積、局所／大域I/F、状態と回復を仕様v1.0へ統合した。R4の成立根拠に加え、後退・テーパ領域の4配置12右辺と4物理荷重問題を照合。構造1.1の採用判断を終了し、現在地は1.2材料・積層・等価入力の統合。構造審査全体・本体実装開始は未完了。'
    fronts=['README.md','rd/planning/preimplementation-plan.md','rd/planning/structural-model-adoption-review.md','rd/planning/requirements-evidence-map.md']
    for name in fronts:
        p=ROOT/name;lines=p.read_text(encoding='utf-8').splitlines();changed=0
        for i,line in enumerate(lines):
            if '参照係数PG-ANCFのR4整合審査を完了し、R5の主シェル採用審査候補として選定。' in line:
                prefix=line.split('参照係数PG-ANCFのR4整合審査を完了し、R5の主シェル採用審査候補として選定。')[0]
                base='docs/' if name=='README.md' else ''
                lines[i]=prefix+common+f' [R5採用判定]({base}verification/paired-shell-r5-results.md)、[採用仕様]({base}main-shell-adopted-specification.md)、[現在地]({base}roadmap.md)。';changed+=1
        assert changed==1,(name,changed)
        if name=='README.md':
            lines=[('**現在地：実装前6作業の「1. 構造モデル採用」→1.2 材料・積層・等価入力の統合。** 主シェルは採用済み。構造審査全体は未完了。C++23の両OSビルド・CLI試験の基盤はあるが、製品の構造・空力・連成計算本体は未実装。' if s.startswith('**現在地：') else s) for s in lines]
        p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'docs/specs/design-decisions.md';lines=p.read_text(encoding='utf-8').splitlines()
    for i,line in enumerate(lines):
        if line.startswith('| ST-01 |') or line.startswith('| ST-04 |'):
            cells=line.split('|')
            content=('P2・R5完了、主シェルに参照係数PG-ANCFを採用。SH-D01〜08、固定full V・Cartesian係数・参照グラフ／領域・必要微分・全トレース・適用／拒否条件を仕様v1.0として固定。'
                if line.startswith('| ST-01 |') else 'P2・R5完了。全3Dエネルギー、係数対のP7補完、参照だけの条件付け、全列24-jet QR、層／領域／背景セル交点で分ける求積を採用。ANS/EAS・法線Schur・人工安定化は既定に追加しない。')
            cells[-2]=' '+content+' 主シェルの採用判断は終了。構造全体・本体開始は未完了（[採用仕様](main-shell-adopted-specification.md)、[R5判定](verification/paired-shell-r5-results.md)、[後続](structural-integration-after-shell.md)） ';lines[i]='|'.join(cells)
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'rd/planning/preimplementation-history.md';lines=p.read_text(encoding='utf-8').splitlines()
    for i,line in enumerate(lines):
        if line.startswith('**現在地は1.1主シェル、R5'):lines[i]='**現在地は1.2材料・積層・等価入力の統合。主シェル1.1／R5は完了。** '+common+' [R5採用判定](verification/paired-shell-r5-results.md)、[採用仕様](main-shell-adopted-specification.md)、[後続の順序・終了条件](structural-integration-after-shell.md)。'
        elif line=='### 主シェルを閉じるまで':lines[i]='### 主シェルの採用まで：R1〜R5完了'
        elif line.startswith('| **R4 '):lines[i]=line.replace('製品全受入・主シェル最終採用・本体開始は未完了','R4単独で主シェル最終採用とはしない。最終判断は次行のR5')
        elif line.startswith('| **R5 主シェル'):lines[i]='| **R5 主シェルの採用判断と仕様固定：完了** | R4の根拠に後退・テーパ領域の品質／物理荷重を追加。適用／拒否、必要微分、局所／大域、状態・接続・回復・費用を[採用仕様v1.0](main-shell-adopted-specification.md)へ統合 | [参照係数PG-ANCFを主シェルに採用](verification/paired-shell-r5-results.md)。1.1を閉じ、1.2→1.3→1.4→1.5へ渡す。構造全体・本体開始は未完了 |'
        elif line.startswith('R1〜R5は実装前6作業'):lines[i]='R1〜R5は実装前6作業の「1」の中で主シェルを閉じる作業で、今回R5まで完了した。次は[1.2材料・積層・等価入力の統合](structural-integration-after-shell.md)。主シェルの候補探索を根拠なく再開しない。項目数を製品完成率や採用日に換算しない。'
        elif line.startswith('- 柱面R2/R3・R4.1と、今回の参照係数R4'):lines[i]='- R4/R5の成立範囲と採用仕様を再利用し、旧PG-Ritzの不採用は維持する。次は1.2の残る材料統合・等価入力を閉じる。新しい具体的反例がない限り、完了済みの主シェル基礎検査へ戻らない。'
        elif line.startswith('| **1 構造モデル採用**'):lines[i]='| **1 構造モデル採用** | **進行中・主シェル1.1採用済み**。R4代表整合・R5仕様固定を終了 | **1.2材料・等価入力→1.3梁→1.4接続・支持・慣性→1.5荷重・回復・構造統合** | ST-01〜10/CV-03の採用式・適用条件・必要微分・接続・独立参照と検証計画が揃う |'
        elif line.startswith('| 1.1 主シェル'):lines[i]='| 1.1 主シェルの補間・ロッキング・求積・幾何 | **R5完了・参照係数PG-ANCFを採用**。[仕様v1.0](main-shell-adopted-specification.md)と[R5採否](verification/paired-shell-r5-results.md) | 主シェルの採用判断は終了。実務入力・全接続・製品性能の担当と条件を後続へ明示して引継ぎ |'
        elif line.startswith('| 1.2 材料・積層'):lines[i]='| **1.2 材料・積層・初期ひずみ・非線形微分：現在地** | 材料点MAT-D04〜07とR4の同じSの積層非線形・層分割・両側回復を再利用 | 全3D主経路との統合、STR-01等価入力の許容形式・不足情報・質量／断面仕事の対応を固定。実材料の適用範囲は受入ケース5へ条件を渡す |'
        elif line.startswith('| 1.4 質量・偏心慣性'):lines[i]='| 1.4 接続・支持・質量・偏心慣性 | 全トレース・厚さ変換等の部分契約と証拠あり。全体は未完了 | 梁―シェル・折れ・交差・厚さ段差・剛体・指定運動、独立拘束・反力・慣性を同じ運動学から閉じる |'
        elif line.startswith('| 1.5 構造の採用統合'):lines[i]='| 1.5 荷重・減衰・結果回復と構造採用統合 | 未完了 | 1.1〜1.4と追従荷重／接線・支持／物理減衰・切断面回復を統合し、ST系の終了条件を満たして2へ渡す |'
        elif line.startswith('1.2〜1.4は1.1の全実験'):lines[i]='1.1の採用判断を閉じたため、共有未知量・材料・接続を同じ主シェルへ統合できる。既存の材料点とR4の証拠は再利用し、次の有限な比較は1.2の未確定事項へ限定する。'
        elif line.startswith('現段階で根拠ある完了率'):lines[i]='主シェルの採用判断は完了した。梁・接続・J1・数値基盤・実務参照と性能の残作業があるため、全体の完了率や日付はまだ根拠を持って算出できない。文書数・試験数を製品完成率へ換算しない。'
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    # Recheck evidence, all links introduced by the R5 status, and final state flags.
    for group in ['sources','artifacts']:
        for name,h in a[group].items():assert sha(ROOT/name)==h,name
    allfronts=fronts+['docs/specs/design-decisions.md','rd/planning/preimplementation-history.md'];links=0
    for name in allfronts:
        p=ROOT/name;s=p.read_text(encoding='utf-8')
        assert '主シェル' in s and ('採用済み' in s or '主シェルとして採用' in s or 'R5完了' in s),name
        for target in re.findall(r'\]\(([^)]+)\)',s):
            if any(key in target for key in ['paired-shell-r5','main-shell-adopted','structural-integration-after-shell']):
                assert (p.parent/target.split('#')[0]).exists(),(name,target);links+=1
    immutable=['docs/specs/structure/main-shell-adopted-specification.md','rd/planning/structural-integration-after-shell.md',
        'rd/structure/main-shell/paired-geometry/records/paired-shell-r5-disposition.json','rd/structure/main-shell/paired-geometry/records/paired-shell-r5-results.md','rd/structure/main-shell/paired-geometry/records/paired-shell-r5-audit.json']
    closure=dict(passed=True,R5_complete=True,main_shell_adopted=True,structural_review_complete=False,product_implementation_start=False,
        current_work='1.2 material/layer/equivalent-input integration',checked_status_links=links,immutable={name:sha(ROOT/name) for name in immutable},
        status_sync_source={Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))},mutable_status_entries=allfronts)
    path=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-shell-r5-closure.json';assert not path.exists();path.write_text(json.dumps(closure,indent=2),encoding='utf-8');print(json.dumps(closure,indent=2))

if __name__=='__main__':main()
