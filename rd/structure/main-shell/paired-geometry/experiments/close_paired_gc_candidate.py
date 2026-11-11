"""Publish the bounded candidate decision only after all qualification gates pass."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
ROOT=rd_paths.ROOT

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    path=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-candidate-result.json';r=json.loads(path.read_text());basis=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-basis-audit.json';b=json.loads(basis.read_text())
    spatial_path=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-spatial-result.json';spatial=json.loads(spatial_path.read_text())
    assert r['adoption_candidate_qualified'] and r['all_numerical_checks'] and b['passed'] and spatial['passed']
    dest=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-candidate-results.md';assert not dest.exists()
    rational=[a for a in r['records'] if a['mode']=='rational'];accuracy=[a for a in r['accuracy'] if a['mode']=='rational']
    emax=max(max(a['independent']['error']) for a in rational);pmax=max(max(a['independent']['position']) for a in rational);ratio=max(a['ratio'] for a in accuracy if a['field']=='error');pratio=max(a['ratio'] for a in accuracy if a['field']=='position')
    rows=[]
    for a in r['records']:
        rows.append(f"| {a['mode']} | {a['amplitude']:g} | {a['n']} | {a['h']:g} | {max(a['independent']['error']):.6g} | {max(a['independent']['position']):.6g} | {a['fit']['rcond']:.4g} | {a['seconds']:.1f} |")
    text=f'''# PG-ANCF：採用審査へ進める候補の成立

2026-11-04。**PG-ANCFを主シェルの最終採用審査へ進める候補として選定する。主シェル最終採用と本体実装開始の判断はまだ行わない。** R2の具体案とR3の今回の決定的な有限厚さ検査を閉じ、次はR4の一般面・物理PDE・実務力学との整合へ進む。

## 何が変わったか

[定式化と主要処理境界](../paired-geometry-ancf-definition.md)の通り、位置と厚さ勾配を別々に補完する代わりに、同じ曲率からφ₁=−κC−zI₀、φ₃=b−κxzI₀を一緒に作る。各スカラー位置形状関数に全3Cartesian係数を掛ける。線形ANCF位置写像、一定の整合質量、変形後の回転に対する閉包、基本r/d/cの独立性と全3D生エネルギーを保持する。実位置が全厚さで連続し、独立な人工面変数の釣合いに依存しない。

この変更で薄肉増幅の未達を変えた。旧GCの独立r/d補完より誤差が大きくなる場合もあるが、1%と2倍を同時に満たすことを条件とし、誤差の小ささだけで採用しない。次数・求積・許容差は追加曲率に合わせて変更していない。[客観性・仕事・固定hの剛体核と柱面の近似上界](../paired-geometry-ancf-invariants.md)も式で整理した。柱面の上界から2倍条件や一般面の性能を結論せず、それぞれ別に判定する。

## 判定した範囲

[事前仕様](paired-geometry-ancf-qualification-protocol.md)の**12配置・計36右辺**を通過。標準3場は正則／歪み、n2/n4、h=.01/.0001の8配置。追加曲率k=1/(1+a²x²)、a=2,4,6は正則／歪みn4、2厚さの4配置。全て解析放物柱面κ=.4、均質等方弾性、片端固定の投影診断である。物理体積力と面力から独立に組み立てたPDE平衡の受入ではない。

追加曲率の最細n4における最大エネルギー相対誤差 **{emax:.9g}**、最大位置相対誤差 **{pmax:.9g}**。薄肉増幅はエネルギー最大 **{ratio:.9g}倍**、位置最大 **{pratio:.9g}倍**。同じ基底で追加曲率も判定し、標準場が数値床に入ったことだけを改善の証拠にしていない。特にa=6の非零誤差と増幅が有効な比較になる。

[補足仕様](paired-geometry-ancf-spatial-protocol.md)に従い、追加曲率のn2・正則／歪み・2厚さの4配置も検査した。**主・補足合計16配置・48右辺**。数値床より大きい追加曲率の誤差でもn2→n4が減少し、両格子の薄肉増幅2倍条件を通過した。[補足結果](paired-geometry-ancf-spatial-result.json)。標準場の数値床以下の増減を収束率の証拠にしない。

| 場 | 歪み振幅 | n | h | 最大エネルギー誤差 | 最大位置誤差 | rcond | 秒 |
|---|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

全数値検査（元点停留、線形残差、全列圧縮、独立15×9観測、支持／共有、質量正値、基準場の元Jacobian照合）を通過。代表自由パッチは剛体6モードのみ。一定F・有限剛体回転・体積と全外面の一定応力仕事・材料＋幾何接線、代表薄肉n2の別SVDも通過。[機械可読の全結果](paired-geometry-ancf-candidate-result.json)。

主系列{r['seconds']:.1f}秒、ピーク{r['peak_bytes']/1024**3:.3f}GiB。事前の1配置300秒・系列1800秒・2GiB以内。補足系列{spatial['seconds']:.1f}秒も別途の300秒・2GiB枠内。これは試作の資源検査であり製品性能の受入ではない。

[基底の独立照合](paired-geometry-ancf-basis-audit.json)では元の原始関数との値差最大{max(x['value'] for x in b['representation']):.3g}、導関数差最大{max(x['derivative'] for x in b['representation']):.3g}。基準曲げ場のS空間への包含は8条件で位置差最大{max(x['position'] for x in b['native_inclusion']):.3g}、ひずみ差最大{max(x['strain'] for x in b['native_inclusion']):.3g}。極小な標準誤差は基準場と補正を分けた評価の値であり、この表現照合限界以下の実変位精度を主張しない。一般近似定理・一様inf-supや、人による独立レビューの実施も主張しない。

## 不採用・停止を保持する

| 仕様 | 判断と根拠 |
|---|---|
| MMS-1 | [有限厚さの薄肉増幅で不採用](mixed-mean-stress-finite-curved-results.md)。[R1原因診断](mixed-mean-stress-finite-cause-results.md)の共同持上げ条件を引き継ぐ |
| C1-WM 2案 | [v1](c1-work-mean-v1-disposition.json)・[v2](c1-work-mean-v2-disposition.json)とも精度未達。延長しない |
| GC v1〜v3 | 連続性と曲率／仕事の閉包を修正して比較したが、薄肉増幅未達。[v1](geometry-compatible-ancf-v1-disposition.json)、[v2](geometry-compatible-ancf-v2-disposition.json)、[v3](geometry-compatible-ancf-v3-disposition.json) |
| GC v4〜v5 | [低求積のaliasと補助rank](geometry-compatible-ancf-v4-screen-disposition.json)、[基底rank停止](geometry-compatible-ancf-v4-basis-stop.json)、[境界準備停止](geometry-compatible-ancf-v5-boundary-stop.json)を保持。合格へ読み替えない |
| GC v6 | [元系列のメモリ超過](geometry-compatible-ancf-candidate-result.json)は同一空間の逐次評価で解消し、[物理場一致](geometry-compatible-ancf-stream-equivalence.json)を確認。ただし[追加曲率の増幅は未達](geometry-compatible-ancf-candidate-stream-result.json)で不採用 |
| セル境界への曲率区間変更 | [2厚さの棄却比較](geometry-compatible-ancf-cell-knots-screen.json)でも増幅未達。全面検証へ進めない |
| PGの対の位置補完 | 上記と異なる独立座標の関係を定義し、同じ3Dエネルギー・基準で今回の適格性検査を通過 |

## 主シェルの最終採用へ残す判断

1. 一般曲面の参照Ritz持上げ、二階微分、チャート接続・全トレース拘束。柱面の式だけで一般面の品質を保証しない。
2. 独立な体積力／面力からの物理PDE求解、反力と仕事、既存の公開参照の適用条件との照合。
3. 採用済み材料点契約と層別求積／回復・初期ひずみを接続し、非線形U/f/K・慣性／固有特性・状態再開を同じ座標で閉じる。
4. 広い支持の補完と疎／境界付き行列の容量・前処理・実務規模の費用を確認する。小規模の資源枠通過から製品性能を推定しない。

R4では式・I/F・不変条件と採用を左右する代表照合を閉じ、網羅的な製品検証は実装と対応させる。これらを一つの仕様へ対応させ、R5で採用／不採用を判断する。要件・閾値を変える必要がある場合を除き、確認待ちを挟む必要はない。[全体ロードマップ](../roadmap.md)、[判断と出典](paired-geometry-ancf-candidate-disposition.json)、[監査](paired-geometry-ancf-candidate-audit.json)。
'''
    dest.write_text(text,encoding='utf-8')
    message='PG-ANCF（位置・厚さ勾配の対のスカラー補完）を採用審査へ進める候補として選定。R2の具体案とR3の主・補足16配置・48右辺の検査を終了。次はR4の一般面・物理PDE・材料／慣性等との整合とR5の最終採用判断。旧MMS-1と独立r/d補完GCの不採用は維持。主シェル最終採用・本体実装開始は未承認。'
    entries=[('README.md','2026-11-04：','docs/'),('rd/planning/preimplementation-plan.md','現在の進行（2026-11-04）：',''),('rd/planning/structural-model-adoption-review.md','現在の到達点（2026-11-04）：',''),('rd/planning/requirements-evidence-map.md','E-SHELLの現在の追加証拠（2026-11-04）：','')]
    for name,prefix,relative in entries:
        p=ROOT/name;lines=p.read_text(encoding='utf-8-sig').splitlines();hits=[i for i,s in enumerate(lines) if s.startswith(prefix)];assert len(hits)==1;lines[hits[0]]=prefix+message+f'[候補判定]({relative}verification/paired-geometry-ancf-candidate-results.md)、[現在地]({relative}roadmap.md)。';p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'docs/specs/design-decisions.md';lines=p.read_text(encoding='utf-8-sig').splitlines()
    for i,line in enumerate(lines):
        if line.startswith('| ST-01 |') or line.startswith('| ST-04 |'):
            cols=line.split('|');cols[-2]=' P2・'+message+'（[定式化](paired-geometry-ancf-definition.md)、[候補判定](verification/paired-geometry-ancf-candidate-results.md)、[計画](roadmap.md)） ';lines[i]='|'.join(cols)
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    # The roadmap is maintained separately so its full domain overview is retained.
    sources=[path,basis,spatial_path,dest,ROOT/'rd/structure/main-shell/paired-geometry/paired-geometry-ancf-definition.md',ROOT/'rd/structure/main-shell/geometry-compatible/geometry-compatible-ancf-definition.md',ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-preparation.json',ROOT/'rd/structure/main-shell/mixed-mean-stress/records/mixed-mean-stress-finite-cause-audit.json',Path(__file__)]
    sources+=list((rd_paths.RECORDS).glob('geometry-compatible-ancf-*.json'))
    sources+=list((rd_paths.RECORDS).glob('c1-work-mean-*-disposition.json'))
    sources+=[ROOT/'rd/structure/main-shell/geometry-compatible/records/geometry-compatible-ancf-cell-knots-protocol.md',ROOT/'rd/structure/main-shell/paired-geometry/experiments/update_paired_gc_roadmap.py',ROOT/'rd/structure/main-shell/paired-geometry/paired-geometry-ancf-invariants.md']
    disposition=ROOT/'rd/structure/main-shell/paired-geometry/records/paired-geometry-ancf-candidate-disposition.json';assert not disposition.exists()
    disposition.write_text(json.dumps(dict(adoption_candidate='PG-ANCF',adoption_candidate_qualified=True,R2_candidate_definition_complete=True,R3_bounded_qualification_complete=True,main_shell_adopted=False,production_ready=False,physical_PDE_acceptance=False,human_independent_review=False,next_stage='R4 general geometry and practical mechanics; R5 main shell adoption',sources={p.relative_to(ROOT).as_posix():digest(p) for p in sources}),indent=2),encoding='utf-8')
    print('Candidate report, decision and current entrypoints updated.')

if __name__=='__main__':main()
