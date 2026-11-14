"""Reconcile mutable status after material 1.2 adoption; keep evidence frozen."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,re
ROOT=rd_paths.ROOT
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    path=ROOT/'rd/structure/material/records/material-shell-adoption-audit.json';audit=json.loads(path.read_text());assert audit['passed'] and audit['material_1_2_complete']
    for group in ['sources','artifacts']:
        for name,h in audit[group].items():assert sha(ROOT/name)==h,name
    common='1.2材料・積層・等価入力の数値仕様と主シェル結合を採用。等方・3D直交異方性・明示全3D有効則を同じ層別U/f/K/Mへ接続し、初期ひずみ・参照軸・全厚・両側回復・密度／慣性を統合した。追加10群と28入力負対照が成立。主シェル1.1／R5と1.2は完了し、現在地は1.3梁の採用。実材料の用途別受入は5へ条件を明示して引継ぎ、構造全体・本体実装開始は未完了。'
    fronts=['README.md','rd/planning/preimplementation-plan.md','rd/planning/structural-model-adoption-review.md','rd/planning/requirements-evidence-map.md']
    for name in fronts:
        p=ROOT/name;lines=p.read_text(encoding='utf-8').splitlines();count=0
        for i,line in enumerate(lines):
            if 'R5を完了し、参照係数PG-ANCFを主シェルとして採用。' in line:
                prefix=line.split('R5を完了し、参照係数PG-ANCFを主シェルとして採用。')[0];base='docs/' if name=='README.md' else ''
                lines[i]=prefix+common+f' [1.2採用判定]({base}verification/material-shell-adoption-results.md)、[材料採用仕様]({base}material-shell-adopted-specification.md)、[現在地]({base}roadmap.md)。';count+=1
            if name=='README.md' and line.startswith('**現在地：'):
                lines[i]='**現在地：実装前6作業の「1. 構造モデル採用」→1.3 梁の採用。** 主シェル1.1と材料1.2は採用済み。構造審査全体は未完了。C++23の両OSビルド・CLI試験の基盤はあるが、製品の構造・空力・連成計算本体は未実装。'
        assert count==1,(name,count);p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'docs/specs/design-decisions.md';lines=p.read_text(encoding='utf-8').splitlines()
    updates={
        'CV-03':'P3・主シェルの材料軸・層順・Mandel／工学変換・Green／第二PK・界面両側・出力タグをMAT-D08〜13として採用。10群照合と28負対照を完了。梁と出力ファイルの具体型は1.3／4で統合',
        'ST-03':'P2・1.2の数値仕様・主シェル結合を採用。全3D等方・直交異方性・明示SPD有効則、層別C/E₀/ρ/参照軸、等価入力と断面仕事・質量の対応を固定。主経路に法線Schurを入れない。実材料の用途別受入は5',
        'ST-05':'P2・採用PG-ANCFの一定整合質量と層別密度を採用。R4と1.2で独立物理Xによる空間慣性、全3D有効入力・断面質量を照合。梁・集中質量・偏心慣性の統合は1.3／1.4に残る',
        'ST-10':'P2/P3・採用PG-ANCFと全3D材料の同じWからU/f/K・両側回復を得る式と処理境界を採用。材料＋幾何接線・独立テンソル経路・方向微分を照合。梁・断面力等の全体回復は1.3〜1.5、生成方式・具体型は4で閉じる'}
    for i,line in enumerate(lines):
        for key,text in updates.items():
            if line.startswith('| '+key+' |'):
                cells=line.split('|');cells[-2]=' '+text+'（[材料採用仕様](material-shell-adopted-specification.md)、[1.2判定](verification/material-shell-adoption-results.md)） ';lines[i]='|'.join(cells)
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'rd/planning/preimplementation-history.md';lines=p.read_text(encoding='utf-8').splitlines()
    for i,line in enumerate(lines):
        if line.startswith('**現在地は1.2材料'):lines[i]='**現在地は1.3梁の採用。主シェル1.1／R5、材料1.2は完了。** '+common+' [1.2採用判定](verification/material-shell-adoption-results.md)、[材料採用仕様](material-shell-adopted-specification.md)、[主シェル採用仕様](main-shell-adopted-specification.md)。'
        elif line.startswith('R1〜R5は実装前6作業'):lines[i]='R1〜R5は実装前6作業の「1」の中で主シェルを閉じる作業で、R5まで完了した。続く[1.2材料・積層・等価入力](verification/material-shell-adoption-results.md)も採用済み。次は1.3梁の採用。主シェルの候補探索を根拠なく再開しない。項目数を製品完成率や採用日に換算しない。'
        elif line.startswith('- R4/R5の成立範囲と採用仕様'):lines[i]='- R4/R5と材料1.2の成立範囲・採用仕様を再利用する。次は1.3梁の理論・DOF・断面・エネルギー・質量・シェル受渡しを一組として選定する。新しい具体的反例がない限り、完了済みの主シェル基礎検査へ戻らない。'
        elif line.startswith('| **1 構造モデル採用**'):lines[i]='| **1 構造モデル採用** | **進行中・主シェル1.1と材料1.2を採用済み** | **1.3梁→1.4接続・支持・慣性→1.5荷重・回復・構造統合** | ST-01〜10/CV-03の採用式・適用条件・必要微分・接続・独立参照と検証計画が揃う |'
        elif line.startswith('| **1.2 材料・積層'):lines[i]='| 1.2 材料・積層・初期ひずみ・非線形微分 | **採用・統合完了**。MAT-D08〜13、等価3D入力、断面仕事／質量、変厚・層別U/f/K/M、両側回復を[仕様v1.0](material-shell-adopted-specification.md)へ固定 | 実材料・用途条件・有効範囲・参照不確かさは5に明示して引継ぎ。法線Schurを既定経路へ入れる作業は要求しない |'
        elif line.startswith('| 1.3 梁・シェル接続'):lines[i]='| **1.3 梁・シェル受渡し：現在地** | 具体方式は未選定。主シェル・材料の境界は固定済み | 軸／せん断／ねじり・大回転・断面／反り、DOF、内力・全接線・質量・独立参照を一組として選定し、シェルへ渡す位置／方向・仕事を閉じる |'
        elif line.startswith('1.1の採用判断を閉じたため'):lines[i]='1.1と1.2の採用判断を閉じたため、同じ主シェルと材料に梁・接続を統合できる。次の有限な比較は1.3の梁候補と断面・運動学の決定へ限定し、既存の主シェル・材料試験は再利用する。'
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    # Verify immutable evidence and all newly introduced current-entry links.
    for group in ['sources','artifacts']:
        for name,h in audit[group].items():assert sha(ROOT/name)==h,name
    allfronts=fronts+['docs/specs/design-decisions.md','rd/planning/preimplementation-history.md'];links=0
    for name in allfronts:
        p=ROOT/name;s=p.read_text(encoding='utf-8');assert '1.3' in s,name
        for target in re.findall(r'\]\(([^)]+)\)',s):
            if 'material-shell' in target:
                assert (p.parent/target.split('#')[0]).exists(),(name,target);links+=1
    frozen=['docs/specs/structure/material-shell-adopted-specification.md','rd/structure/material/records/material-shell-adoption-disposition.json','rd/structure/material/records/material-shell-adoption-results.md','rd/structure/material/records/material-shell-adoption-audit.json']
    closure=dict(passed=True,material_1_2_complete=True,main_shell_adopted=True,structural_review_complete=False,product_implementation_start=False,
        current_work='1.3 beam adoption and shell interface',immutable={p:sha(ROOT/p) for p in frozen},
        status_sync_source={Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))},mutable_status_entries=allfronts,checked_status_links=links)
    dest=ROOT/'rd/structure/material/records/material-shell-adoption-closure.json';assert not dest.exists();dest.write_text(json.dumps(closure,indent=2),encoding='utf-8');print(json.dumps(closure,indent=2))

if __name__=='__main__':main()
