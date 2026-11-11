"""Publish a successful R4 evidence disposition and reconcile six status entries."""
from pathlib import Path
import hashlib,json,re
ROOT=Path(__file__).resolve().parents[2]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    audit_path=ROOT/'docs/verification/paired-reference-r4-audit.json';audit=json.loads(audit_path.read_text());assert audit['passed']
    for section in ['sources','artifacts']:
        for path,digest in audit[section].items():assert sha(ROOT/path)==digest,path
    m=audit['metrics'];dest=ROOT/'docs/verification/paired-reference-r4-disposition.json';assert not dest.exists()
    disposition=dict(R4_complete=True,selected_for_R5='reference-coefficient PG-ANCF with exact integral paired lifting',
        rejected_membrane_Ritz_binding_remains_rejected=True,main_shell_adopted=False,structural_review_complete=False,
        product_implementation_start=False,metrics=m,remaining_milestone='R5: integrate applicability/rejection conditions and freeze main-shell adoption specification',
        limitations=['representative R4 evidence, not all practical application acceptance',
            '100800-coordinate capacity measured on independent supported two-curved-patch blocks; globally connected product performance not measured',
            'exact modal restart oracle is not product time integrator adoption',
            'prototype coefficient adapter is not product parallel runtime'],
        sources={'docs/verification/paired-reference-r4-audit.json':sha(audit_path),Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))})
    dest.write_text(json.dumps(disposition,indent=2),encoding='utf-8')
    report=ROOT/'docs/verification/paired-reference-r4-results.md'
    report.write_text(f'''# 参照係数PG-ANCF：R4完了と採用審査候補の選定

2026-11-05。**R4「実務に必要な力学との整合」を完了し、参照係数PG-ANCFをR5の主シェル採用審査候補として選定する。** 主シェルの最終採用、構造審査全体、製品計算本体の実装開始は未完了。次はR5の適用／拒否条件の統合と仕様固定である。[採否JSON](paired-reference-r4-disposition.json)、[独立集約監査](paired-reference-r4-audit.json)。

## 何を変えて、問題を解消したか

旧[PG-Ritz膜最小化への固定案](pg-ritz-r4-results.md)の不採用は維持する。今回の候補は参照Ritz解に物理空間を固定せず、曲面の多項式係数から位置・厚さ勾配の対を直接積分する。[定式化](../paired-reference-coefficient-definition.md)。膜の残る成分も全3Dエネルギーに含める。これはRitz持上げの成功を主張するものではなく、その置換案の採用審査である。

有理係数で空間・必須場・全トレースを構成し、基本directorをQ2/√gとする。参照の離散幾何内積による一定の座標変換を使い、係数の組合せを高精度で完了してから倍精度へ変換する。連続積分全体の区間保証と、求積行列の合同変換の保証は区別する。位置や剛性の小さい特異方向を捨てていない。現在配置や材料・荷重に依存しないNと一定参照質量を維持する。

途中の[独立成分追加案](coefficient-compatible-v1-entrance-result.json)、[係数対の初期版](exact-paired-primary-v2-result.json)、[区間Choleskyの旧実装](mass-paired-v1-precision-result.json)は上書きせず保持した。初期版の未達を、現在の候補の合格で取り消さない。

## 同じ候補で揃えた証拠

| 審査群 | 実施結果 |
|---|---|
| 一般曲面の定式化 | [有理数の恒等式](paired-lift-exact-result.json)：全4次単項式を含む曲面、幅・歪み2チャート、両方向256場。膜2条件、横せん断／法線条件、二階混合微分が係数上で厳密に零。残る膜成分は消さない |
| 決定的有限厚さ | [主・非二次品質](mass-paired-quality-result.json)24配置72右辺と[追加係数入力](mass-paired-coefficient-input-result.json)4配置12右辺、計28配置84右辺。柱面・楕円・ねじれ・鞍面、真の非二次・混合4次面とアフィン歪み、2厚さ・2解像度。最細1%・増幅2倍・細分化の全条件を通過 |
| 数値と参照 | [参照検査](mass-paired-reference-result.json)：真の法線厚さ、全域の正則性十分条件、独立複素微分、参照再現、全トレース核。全品質配置で剛性／質量平方根の3閾値rank、原停留、条件数、固定端、全3D作用素の同値性を通過 |
| 独立な物理荷重 | [非アフィン体積力・全面力](mass-paired-body-result.json)：u_x=aX²から独立にb=−divσとt=σnを定義。非二次面・n1/2・2厚さの4問題で物理誤差、原平衡、外力仕事、全面力釣合いが成立。Kqを外力の生成に使っていない |
| 実際のパッチ接続 | [二曲面パッチの物理PDE](mass-paired-interface-result.json)：同分割1/1と非整合1/2。全トレース一致、界面の転置仕事、独立面力からの平衡解を確認。平均トレースだけの接続ではない |
| 材料・非線形・慣性・状態 | [同じfull Sによる6群](mass-paired-mechanics-result.json)：異方性2層・初期ひずみ／応力、U/f/K、客観性、層分割、両側回復、Newton復帰、整合質量、6剛体、物理Xによる慣性、固有対・負固有値保持、全modeの確定／棄却／再開と16負対照が成立 |
| 求積と向き | [独立求積・法線／層反転](mass-paired-orientation-result.json)：非アフィンU/f/K/Mの45→49比較、層順の同時反転、全6面の仮想仕事。層順の誤りと測度の符号誤りを負対照で検出 |
| 局所／大域と容量 | [主要I/F・原KKT](../paired-reference-local-global-contract.md)を定義。[実容量](mass-paired-capacity-result.json)：100800座標の因子を実際に保持し、作用・求解を確認。独立部材集合の実験であり、一体接続した製品性能の合格ではない |

主・補足品質での最細最大相対エネルギー誤差は {m['finest_energy_error']:.6g}、位置誤差は {m['finest_position_error']:.6g}。薄肉増幅の最大は {m['max_thin_amplification']:.6g} 倍で、2倍以内。原停留の最大は {m['max_original_stationarity']:.6g}。参照準備を含む主・補足品質の最長1配置は {m['max_cold_quality_case_seconds']:.3f} 秒で、300秒以内。数値床1e−12を緩めていない。

容量試験のピークRSSは {m['capacity_peak_rss']/1024**3:.3f} GiB。全結合10万自由度の疎因子fill・前処理反復数・製品フラッタ探索時間は、この数値から推定して合格扱いしない。

[集約条件](paired-reference-r4-protocol.md)に対する独立監査で、9結果群・{audit['counts']['sources']}ソース／参照・{audit['counts']['artifacts']}保存場のハッシュと受入算術を検査した。個々の結果のR4_complete=falseは、その単独結果では全体の完了を宣言しないための記録である。今回の全体判断は本書と採否JSONが行う。

## R5へ渡すもの

1. 参照係数空間、必要微分、実厚さ、全トレースと物理面仕事、局所／大域の処理契約。
2. 28配置84右辺の品質と、独立物理荷重、材料・慣性・状態、実容量の再現可能な保存場。
3. 適用条件と拒否条件。多項式係数を持つ正則な曲面チャートと正の法線厚さ管が前提。任意の点別C² jetだけで準備できるとはしない。幾何変換・近似誤差、厚さ段差・交差接続、外部格子との対応は各統合契約で扱う。

R5で主シェルの採否と仕様を固定し、その後に構造1.2〜1.5を統合する。製品用疎解法・時間積分・強連成・全実務受入の実装結果までを、R4完了と呼ばない。一方、今回不合格だった基礎課題を未解決のままR5へ送ったわけでもない。[現在地と全体順序](../roadmap.md)。
''',encoding='utf-8')
    common='参照係数PG-ANCFのR4整合審査を完了し、R5の主シェル採用審査候補として選定。主有限厚さ16配置と一般非二次12配置、計28配置84右辺で1%・薄肉増幅2倍・rank・原停留が成立。独立物理荷重、非整合二曲面接続、積層非線形・質量・状態・実容量も照合した。旧膜Ritz固定案の不採用は維持。現在地はR5の適用／拒否条件の統合と仕様固定で、主シェル最終採用・構造審査全体・本体実装開始は未完了。'
    fronts=['README.md','docs/preimplementation-plan.md','docs/structural-model-adoption-review.md','docs/requirements-evidence-map.md']
    for name in fronts:
        path=ROOT/name;s=path.read_text(encoding='utf-8');lines=s.splitlines();changed=0
        for i,line in enumerate(lines):
            if 'PG-Ritz一般曲面案を検証し、不採用と判定。' in line:
                prefix=line.split('PG-Ritz一般曲面案を検証し、不採用と判定。')[0]
                link='docs/' if name=='README.md' else ''
                lines[i]=prefix+common+f' [R4完了判定]({link}verification/paired-reference-r4-results.md)、[現在地]({link}roadmap.md)。';changed+=1
        assert changed==1,(name,changed)
        s='\n'.join(lines)+'\n'
        if name=='README.md':s=s.replace('→主シェル選定。','→主シェルR5採用審査・仕様固定。')
        path.write_text(s,encoding='utf-8')
    path=ROOT/'docs/design-decisions.md';lines=path.read_text(encoding='utf-8').splitlines()
    for i,line in enumerate(lines):
        if line.startswith('| ST-01 |') or line.startswith('| ST-04 |'):
            cells=line.split('|');cells[-2]=' P2・参照係数PG-ANCFのR4代表整合審査を完了。28配置84右辺、全トレース・物理荷重・材料・慣性・状態・実容量が成立。R5で主シェル採否と仕様を固定する。旧膜Ritz固定案は不採用のまま。ST全体の基準化・本体実装開始は未完了（[R4完了判定](verification/paired-reference-r4-results.md)、[定式化](paired-reference-coefficient-definition.md)、[処理契約](paired-reference-local-global-contract.md)） ';lines[i]='|'.join(cells)
    path.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    path=ROOT/'docs/roadmap.md';lines=path.read_text(encoding='utf-8').splitlines()
    for i,line in enumerate(lines):
        if line.startswith('**現在地は1.1主シェル、R4'):lines[i]='**現在地は1.1主シェル、R5の採用判断・仕様固定。R4は完了。** '+common+' [R4完了判定](verification/paired-reference-r4-results.md)、[局所／大域・主要I/F](paired-reference-local-global-contract.md)。'
        elif line.startswith('| **R4 '):lines[i]='| **R4 実務に必要な力学との整合：完了** | 参照係数PG-ANCFの定式化と二階微分・実厚さ・全トレース、28配置84右辺の決定的品質、独立体積力／全面力、実際の非整合二曲面接続、積層非線形・慣性・状態、局所／大域と実容量を照合。[全体判断](verification/paired-reference-r4-results.md) | 旧膜Ritz固定案は不採用を維持。参照係数候補をR5へ選定。製品全受入・主シェル最終採用・本体開始は未完了 |'
        elif line.startswith('| R5 主シェル'):lines[i]=line.replace('| R5 主シェルの採用判断と仕様固定 |','| **R5 主シェルの採用判断と仕様固定：現在地** |')
        elif line.startswith('R1〜R5は実装前6作業'):lines[i]='R1〜R5は実装前6作業の「1」の中で主シェルを閉じる作業。R4までの代表整合審査を完了した。次のR5では今回の成立根拠と適用／拒否条件を統合し、主シェルの採否・仕様を固定する。項目数を製品完成率や採用日に換算しない。'
        elif line.startswith('- 柱面R2/R3とR4.1は再実行しない。'):lines[i]='- 柱面R2/R3・R4.1と、今回の参照係数R4の成立範囲を再利用する。旧PG-Ritz不採用を維持し、次はR5の採用判断・仕様固定を行う。具体的な新しい反例がない限り、完了済みの基礎検査へ戻らない。'
    path.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    # Final links and frozen evidence are checked after the six current entries change.
    docnames=['docs/verification/paired-reference-r4-results.md','docs/paired-reference-coefficient-definition.md','docs/paired-reference-local-global-contract.md','docs/verification/paired-reference-r4-protocol.md']
    links=0
    for name in docnames:
        path=ROOT/name
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text(encoding='utf-8')):
            if '://' in target or target.startswith('#'):continue
            target=target.split('#')[0];assert (path.parent/target).exists(),(name,target);links+=1
    for section in ['sources','artifacts']:
        for path,digest in audit[section].items():assert sha(ROOT/path)==digest,path
    closure=dict(passed=True,R4_complete=True,main_shell_adopted=False,product_implementation_start=False,checked_new_links=links,
        immutable={name:sha(ROOT/name) for name in docnames+[dest.relative_to(ROOT).as_posix(),audit_path.relative_to(ROOT).as_posix()]},
        mutable_status_entries=['README.md','docs/preimplementation-plan.md','docs/structural-model-adoption-review.md','docs/requirements-evidence-map.md','docs/design-decisions.md','docs/roadmap.md'])
    (ROOT/'docs/verification/paired-reference-r4-closure.json').write_text(json.dumps(closure,indent=2),encoding='utf-8')
    print(json.dumps(closure,indent=2))

if __name__=='__main__':main()
