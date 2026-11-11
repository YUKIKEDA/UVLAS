"""Advance the six live status entrypoints after immutable 1.4 closure."""
from pathlib import Path
import json,hashlib,re
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    auditpath=ROOT/'docs/verification/connection-adoption-audit.json';audit=json.loads(auditpath.read_text());assert audit['passed'] and audit['connection_1_4_complete'];dest=ROOT/'docs/verification/connection-adoption-closure.json';assert not dest.exists()
    def verify():
        for group in ['sources','artifacts']:
            for name,h in audit[group].items():assert sha(ROOT/name)==h,name
    verify();fronts=['README.md','docs/preimplementation-plan.md','docs/structural-model-adoption-review.md','docs/requirements-evidence-map.md','docs/design-decisions.md','docs/roadmap.md'];common='1.4を完了し、全位置接着・明示した剛体取付・指定支持・集中／偏心慣性の数値仕様を採用。梁―シェル実平衡、組合せ固有特性、折れ・T・交差・厚さ段差、独立連続体との非アフィン比較が成立。接合範囲のパッチ分割と厚さ方向の収束条件を固定。主シェル1.1・材料1.2・梁1.3・接続1.4は採用済み。現在地は1.5荷重・減衰・回復・構造統合。構造全体・本体実装開始は未完了。'
    for name in fronts[:4]:
        p=ROOT/name;lines=p.read_text(encoding='utf-8').splitlines();count=0;base='docs/' if name=='README.md' else ''
        for i,line in enumerate(lines):
            if '1.3を完了し、SE(3)' in line:lines[i]=line.split('1.3を完了し、SE(3)')[0]+common+f' [接続1.4の判定]({base}verification/connection-adoption-results.md)、[接続仕様]({base}connection-adopted-specification.md)、[現在地]({base}roadmap.md)。';count+=1
            if name=='README.md' and line.startswith('**現在地：'):lines[i]='**現在地：実装前6作業の「1. 構造モデル採用」→1.5 荷重・減衰・回復・構造統合。** 主シェル1.1・材料1.2・梁1.3・接続1.4は採用済み。構造審査全体は未完了。C++23の両OSビルド・CLI試験の基盤はあるが、製品の構造・空力・連成計算本体は未実装。'
            if name=='README.md' and line.startswith('**計算本体の実装開始条件は未達。**'):lines[i]='**計算本体の実装開始条件は未達。** 構造1.1〜1.4の定式化・主要契約を採用済み。荷重・回復を含む構造全体、UVLM、連成・時間・安定性、数値基盤と主要I/F、実務受入計画を閉じる必要がある。[実装前設計計画](docs/preimplementation-plan.md) に従って各作業の開始条件を統合する。'
        assert count==1,(name,count);p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'docs/design-decisions.md';lines=p.read_text(encoding='utf-8').splitlines();states={
        'ST-05':'P2・シェル一定M・梁M(q)に、粒子／偏心剛体の物理運動エネルギー、M微分と慣性、重力を統合して1.4完了。実接続の整合M/Kと独立固有参照、厚さ分割の質量加法性を確認。動的方式と具体精度型は3／4',
        'ST-07':'P2・1.4完了。点と接線フレーム、重心偏心・重心慣性、M/M微分・慣性力・重力／接線、粒子との区別、MassId／出自／所有領域・追加／置換と二重計上拒否をCN-D05/06として固定',
        'ST-08':'P2・1.4完了。全位置の実接着面、明示的剛体取付、指定支持C/J/Hと混合時間微分、線形拘束／全SE(3)閉路の大域従属、原反力・仕事を採用。実平衡・モード・4接合と非アフィン比較が成立。接合の面内／厚さ収束条件を保持。時間積分の拘束方式は3',
        'ST-10':'P2/P3・主シェル／材料／梁と接続の同じエネルギー・拘束から全接線と反力を導出。非線形KKT／独立消去、乗数H項を照合。荷重・減衰・応力／断面力／支持を含む全体回復は現在1.5、高階慣性微分と具体型は3／4'}
    for i,line in enumerate(lines):
        for key,state in states.items():
            if line.startswith('| '+key+' |'):
                cells=line.split('|');cells[-2]=' '+state+'（[接続仕様](connection-adopted-specification.md)、[1.4判定](verification/connection-adoption-results.md)） ';lines[i]='|'.join(cells)
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'docs/structural-model-adoption-review.md';lines=p.read_text(encoding='utf-8').splitlines()
    for i,line in enumerate(lines):
        if line.startswith('| ST-05/07 '):lines[i]='| ST-05/07 質量・偏心慣性 | シェル一定M・梁M(q)、[接続／偏心慣性の採用](connection-adopted-specification.md)、厚さ分割の質量保存 | **1.4の数値仕様は完了**。時間残差／具体型・性能は3／4、用途別受入は5 | 同じ運動学のM・必要微分・慣性・重力と台帳を固定済み |'
        elif line.startswith('| ST-08 '):lines[i]='| ST-08 界面・支持・剛体結合 | [接続採用仕様](connection-adopted-specification.md)、[1.4の実平衡・モード・有限面・非アフィン比較](verification/connection-adoption-results.md) | **1.4の採用判断は完了**。任意接合の収束・局所応力受入は5、動的拘束は3 | C/J/H・時間微分、独立化と拒否、原反力、面内／厚さ収束条件を固定済み |'
        elif line.startswith('| ST-10 '):lines[i]='| ST-10 内力・接線・結果回復 | シェル・材料・梁・接続の採用エネルギーと全接線、原反力・仕事 | **現在地1.5**。荷重・減衰と応力／断面力／接合・支持の全体回復を統合 | 同じ状態・参照・測度・符号・版の処理仕様と回復適格性を閉じる |'
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'docs/roadmap.md';lines=p.read_text(encoding='utf-8').splitlines()
    for i,line in enumerate(lines):
        if line.startswith('**現在地は1.4接続'):lines[i]='**現在地は1.5荷重・減衰・回復・構造統合。1.1〜1.4は採用済み。** '+common+' [1.4判定](verification/connection-adoption-results.md)、[接続仕様](connection-adopted-specification.md)。'
        elif line.startswith('R1〜R5で主シェル1.1を閉じ'):lines[i]='R1〜R5で主シェル1.1を閉じ、材料1.2・梁1.3・接続1.4も採用済み。次は1.5で荷重・減衰・回復と構造全体の契約を統合する。項目数を製品完成率や採用日に換算しない。'
        elif line.startswith('- 主シェルR4/R5、材料1.2、梁1.3の証拠'):lines[i]='- 主シェルR4/R5・材料1.2・梁1.3・接続1.4の成立根拠を再利用する。1.5は荷重・減衰・回復と構造全体の整合を対象にし、既に閉じた候補比較を再開しない。'
        elif line.startswith('| **1 構造モデル採用**'):lines[i]='| **1 構造モデル採用** | **進行中・1.1主シェル／1.2材料／1.3梁／1.4接続を採用済み** | **現在1.5：荷重・減衰・回復・構造統合** | ST-01〜10/CV-03の採用式・適用条件・必要微分・接続・独立参照と検証計画が揃う |'
        elif line.startswith('| **1.4 接続・支持・質量'):lines[i]='| 1.4 接続・支持・質量・偏心慣性 | **採用完了**。[接続仕様](connection-adopted-specification.md)、[実平衡・モード・有限厚さ接合と慣性の検証](verification/connection-adoption-results.md) | 接合範囲のパッチ分割と面内／厚さ収束を適用条件に保持。動的拘束と具体数値基盤は3／4、実務接合の受入は5 |'
        elif line.startswith('| 1.5 ') or line.startswith('| **1.5 '):lines[i]='| **1.5 荷重・減衰・回復・構造統合：現在地** | 1.1〜1.4のエネルギー・全接線・点／層／断面と反力の回復境界を採用済み。[有限な審査項目](verification/structural-final-integration-protocol.md)を固定 | dead／follower荷重・物理面測度・荷重接線、減衰と散逸、応力／断面力／接合／支持の全体回復、構造の主要処理・状態版・必要微分を統合して1を閉じる |'
        elif line.startswith('1.1〜1.3の採用判断を閉じたため'):lines[i]='1.1〜1.4の採用判断を閉じた。1.5は同じ状態のエネルギー・荷重・減衰・回復を統合し、構造全体の採用式と主要契約を確定する。構造採用を閉じた後に2 UVLMへ進む。'
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8');verify();links=0
    for name in fronts:
        p=ROOT/name;s=p.read_text(encoding='utf-8');assert '1.5' in s and 'connection-adopt' in s,name
        for target in re.findall(r'\]\(([^)]+)\)',s):
            if 'connection-adopt' in target:assert (p.parent/target).exists(),(name,target);links+=1
    frozen=['docs/connection-adopted-specification.md','docs/verification/connection-adoption-disposition.json','docs/verification/connection-adoption-results.md','docs/verification/connection-adoption-audit.json'];closure=dict(passed=True,connection_1_4_complete=True,beam_1_3_complete=True,material_1_2_complete=True,main_shell_adopted=True,structural_review_complete=False,product_implementation_start=False,current_work='1.5 loads damping recovery and structural integration',immutable={name:sha(ROOT/name) for name in frozen},status_sync_source={Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))},mutable_status_entries=fronts,checked_status_links=links);dest.write_text(json.dumps(closure,indent=2),encoding='utf-8');print(json.dumps(closure,indent=2))
if __name__=='__main__':main()
