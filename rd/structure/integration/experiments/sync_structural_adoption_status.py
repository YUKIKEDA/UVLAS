"""Close live structural status and advance current work to UVLM adoption."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,re
ROOT=rd_paths.ROOT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    auditpath=ROOT/'rd/structure/integration/records/structural-adoption-audit.json';a=json.loads(auditpath.read_text());assert a['passed'] and a['structural_review_complete'];dest=ROOT/'rd/structure/integration/records/structural-adoption-closure.json';assert not dest.exists()
    def verify():
        for group in ['sources','artifacts']:
            for name,h in a[group].items():assert sha(ROOT/name)==h,name
    verify();fronts=['README.md','rd/planning/preimplementation-plan.md','rd/planning/structural-model-adoption-review.md','rd/planning/requirements-evidence-map.md','docs/specs/design-decisions.md','rd/planning/preimplementation-history.md'];common='1.5の荷重・減衰・支持・応力／断面力回復・時系列・状態規約を統合し、構造モデル採用「1」を完了。主シェル1.1・材料1.2・梁1.3・接続1.4を維持し、新規25群と実組合せ平衡・全残差接線・独立参照・保存復元を検証。現在地は2 UVLM採用。時間・連成・高階微分は3、具体型／疎解法・性能は4、実務受入は5、実装開始審査は6に残り、本体実装開始条件は未達。'
    for name in fronts[:4]:
        p=ROOT/name;lines=p.read_text(encoding='utf-8').splitlines();count=0;base='docs/' if name=='README.md' else ''
        for i,line in enumerate(lines):
            if '1.4を完了し、全位置接着' in line:lines[i]=line.split('1.4を完了し、全位置接着')[0]+common+f' [構造採用の終了判断]({base}verification/structural-adoption-results.md)、[統合仕様]({base}structural-adopted-specification.md)、[現在地]({base}roadmap.md)。';count+=1
            if name=='README.md' and line.startswith('**現在地：'):lines[i]='**現在地：実装前6作業の「2. UVLM採用」。** 「1. 構造モデル採用」は1.1〜1.5を完了。構造・空力・連成の製品計算本体は未実装で、実装開始審査は6。C++23のWindows／Linuxビルド・CLI試験の基盤は確認済み。'
            if name=='README.md' and line.startswith('**計算本体の実装開始条件は未達。**'):lines[i]='**計算本体の実装開始条件は未達。** 構造の採用式・主要処理・適用境界は確定。次はUVLMの採用、連成・時間・安定性、数値基盤と具体I/F、実務受入計画を閉じ、[実装開始条件](rd/planning/preimplementation-plan.md)を審査する。'
            if name=='rd/planning/preimplementation-plan.md' and line.startswith('| 理論・構造 |'):lines[i]='| 理論・構造 | PG-ANCF・全3D材料／積層・SE(3)梁・接続／慣性・荷重／減衰／回復を[統合採用](structural-adopted-specification.md)。1.1〜1.5完了 | 3／4の動的高階微分・具体型／疎解法、5の実務受入、6の全体開始審査 | **構造モデル採用1は完了** |'
        assert count==1,(name,count);p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'rd/planning/structural-model-adoption-review.md';lines=p.read_text(encoding='utf-8').splitlines()
    for i,line in enumerate(lines):
        if line.startswith('2026-10-31。実装前の残作業1'):lines[i]='2026-10-31に開始した構造モデルの採用審査。**2026-11-05に1.1〜1.5を統合して終了。** 現在の採用式・適用条件は[構造統合仕様](structural-adopted-specification.md)、根拠は[終了判断](verification/structural-adoption-results.md)を参照。以下の初期候補・未達の記録は履歴として保持する。'
        if line.startswith('| ST-02 '):lines[i]='| ST-02 参照配置・幾何 | [主シェル仕様](main-shell-adopted-specification.md)、[構造統合仕様](structural-adopted-specification.md)、接合参照面・後退テーパ | **構造モデル採用は完了**。外部格子の具体変換／誤差予算は4／5 | 参照写像・物理厚さ／法線・正則性・面対応を固定済み |'
        elif line.startswith('| ST-05/07 '):lines[i]='| ST-05/07 質量・偏心慣性 | [接続仕様](connection-adopted-specification.md)、[同一平衡のM/K/C](verification/structural-combined-result.json)、質量台帳 | **構造モデル採用は完了**。動的高階微分・具体型／性能は3／4 | 同じ運動学のM・慣性・重力・所有を統合済み |'
        elif line.startswith('| ST-09 '):lines[i]='| ST-09 荷重・支持・物理減衰 | [統合仕様](structural-adopted-specification.md)、[荷重／減衰](verification/structural-load-damping-result.json)、[時系列と独立支持参照](verification/structural-history-state-result.json) | **1.5採用完了**。時間ステージと全空力依存は2／3、係数の実材料同定は5 | 全荷重接線、客観粘性、動く支持の仕事、イベント／片側／外挿拒否を固定済み |'
        elif line.startswith('| ST-10 '):lines[i]='| ST-10 内力・接線・結果回復 | [統合仕様](structural-adopted-specification.md)、[応力／切断面](verification/structural-recovery-result.json)、[全残差の接線](verification/structural-combined-result.json) | **1.5採用完了**。製品生成方式・具体型は4、用途別回復受入は5 | 同じ状態・測度・軸・層側・切断面・基準点と収支、原反力を固定済み |'
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'docs/specs/design-decisions.md';lines=p.read_text(encoding='utf-8').splitlines();states={
        'ST-01':'P2/P3・SH-D01〜08の参照係数PG-ANCFを採用し構造1を完了。一般曲面・薄肉・物理PDE・後退テーパ・全位置トレースの成立根拠を保持。具体疎解法／精度は4、実務全形状は5',
        'ST-02':'P2/P3・固定参照写像、実法線／厚さ、正則性、後退テーパ領域、接合物理面の対応を採用。構造1完了。外部格子の具体的変換と幾何誤差予算は4／5',
        'ST-03':'P2/P3・MAT採用仕様の全3D Green／第二PK・等価入力・積層・材料軸・初期ひずみ・両側回復を統合。構造1完了。実材料の用途別受入は5',
        'ST-04':'P2/P3・採用PG空間・全3Dエネルギー・求積・参照条件付けの根拠を維持し構造1完了。未達の旧方式や人工補正を既定へ戻さない。具体精度／疎求解は4',
        'ST-05':'P2/P3・シェル一定M、梁／剛体M(q)と慣性、重力、同一平衡M/K/Cを統合し構造1完了。動的接線M,qqと必要三階運動学は3／4へ明示して引継ぎ',
        'ST-06':'P2/P3・BM-D01〜09のSE(3)梁を採用し構造1完了。軸・せん断・ねじり・曲げ、連成断面と初期ひずみ、自由反り・拘束反り拒否・局所回復・精度条件を維持',
        'ST-07':'P2/P3・CN-D05/06の粒子／偏心剛体・M／必要微分・慣性／重力、質量の所有／出自／追加置換と二重計上拒否を統合し構造1完了',
        'ST-08':'P2/P3・CN-D01〜09の全位置接着・明示剛体取付・指定支持・冗長／不整合・反力と仕事を採用し構造1完了。動的拘束と初期化は3、疎構造は4、接合用途別収束は5',
        'ST-09':'P2/P3・SI-D04/05/07として荷重測度・物体／世界成分・圧力／全接線、客観粘性、弾性／減衰支持・基準仕事、明示時系列を採用し1.5完了。時刻／空力依存は2／3、実材料係数は5',
        'ST-10':'P2/P3・SI-D06/08/09として全残差・全接線、測度／軸／層側の応力、物理切断面・原反力・収支、状態規約を統合し1.5完了。必要動的高階微分は3／4、製品生成・型は4、実務回復受入は5',
        'CV-03':'P3・採用材料・梁・構造統合仕様で測度・基底・層／片側・Mandel／工学・断面共役量と物理wrenchを固定。構造1完了。製品出力型とファイル形式は4'}
    for i,line in enumerate(lines):
        for key,value in states.items():
            if line.startswith('| '+key+' |'):
                cells=line.split('|');cells[-2]=' '+value+'（[構造統合仕様](structural-adopted-specification.md)、[終了判断](verification/structural-adoption-results.md)） ';lines[i]='|'.join(cells)
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'rd/planning/preimplementation-history.md';lines=p.read_text(encoding='utf-8').splitlines()
    for i,line in enumerate(lines):
        if line.startswith('**現在地は1.5'):lines[i]='**現在地は2 UVLM採用。「1. 構造モデル採用」は1.1〜1.5を完了。** '+common+' [構造採用の終了判断](verification/structural-adoption-results.md)、[統合仕様](structural-adopted-specification.md)。'
        elif line.startswith('R1〜R5で主シェル1.1を閉じ'):lines[i]='R1〜R5で主シェル1.1を閉じ、材料1.2・梁1.3・接続1.4・荷重／減衰／回復1.5も統合採用した。構造モデル採用1を終了し、次は2 UVLMの方式・精度参照を閉じる。完了項目数を製品完成率へ換算しない。'
        elif line.startswith('- 主シェルR4/R5・材料1.2・梁1.3・接続1.4の成立根拠'):lines[i]='- 構造1.1〜1.5の採用根拠・式・主要契約を再利用する。現在の2 UVLMで構造の候補探索を再開しない。'
        elif line.startswith('    A[要求合意'):lines[i]='    A[要求合意・ビルド基盤] --> B[1 構造採用：完了]'
        elif line.startswith('    B --> C[2 UVLM採用]'):lines[i]='    B --> C[2 UVLM採用：現在地]'
        elif line.startswith('| **1 構造モデル採用**'):lines[i]='| 1 構造モデル採用 | **完了**。1.1〜1.5を[統合採用](structural-adopted-specification.md)、[ST-01〜10/CV-03と新規25群を監査](verification/structural-adoption-results.md) | 時間・高階微分は3、具体数値実装は4、用途別受入は5、全体開始審査は6へ渡す | 構造の採用式・適用／拒否・必要微分・主要契約・独立根拠が対応済み |'
        elif line.startswith('| 2 UVLM採用 |'):lines[i]='| **2 UVLM採用：現在地** | 渦核・荷重・後流時刻・調和応答等の契約と限定証拠あり。精度参照・最終方式は未確定 | 渦配置・コア・自己項・後縁・荷重・後流を一組で選定し、定常／非定常・複数面・幾何・収束を審査 | AE系の採用式・履歴・荷重と独立参照、数値誤差予算・拒否条件が揃う |'
        elif line.startswith('| **1.5 '):lines[i]='| 1.5 荷重・減衰・回復・構造統合 | **採用完了**。[統合仕様](structural-adopted-specification.md)、[終了判断](verification/structural-adoption-results.md)。構造1全体を終了 | 必要な動的高階微分／積分器を3、具体型／疎解法／性能を4、実務受入を5、実装開始を6で審査 |'
        elif line.startswith('1.1〜1.4の採用判断を閉じた。'):lines[i]='1.1〜1.5を閉じ、構造モデル採用1を終了した。現在は2 UVLM。構造の動的・実装・実務受入の未完了は3〜6に明記して継承する。'
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8');verify();links=0
    for name in fronts:
        p=ROOT/name;s=p.read_text(encoding='utf-8');assert 'UVLM' in s and 'structural-adopt' in s,name
        for target in re.findall(r'\]\(([^)]+)\)',s):
            if 'structural-adopt' in target:assert (p.parent/target).exists(),(name,target);links+=1
    frozen=['docs/specs/structure/structural-adopted-specification.md','rd/structure/integration/records/structural-adoption-disposition.json','rd/structure/integration/records/structural-adoption-results.md','rd/structure/integration/records/structural-adoption-audit.json'];closure=dict(passed=True,structural_1_5_complete=True,structural_review_complete=True,product_implementation_start=False,current_work='2 UVLM model adoption',immutable={name:sha(ROOT/name) for name in frozen},status_sync_source={Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))},mutable_status_entries=fronts,checked_status_links=links);dest.write_text(json.dumps(closure,indent=2),encoding='utf-8');print(json.dumps(closure,indent=2))
if __name__=='__main__':main()
