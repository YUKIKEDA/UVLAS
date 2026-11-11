"""Publish bounded R4.1 evidence while retaining the original resource stop."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    path=ROOT/'docs/verification/paired-geometry-ancf-physical-pde-batch-v3-result.json';r=json.loads(path.read_text());assert r['passed'] and len(r['records'])==16
    report=ROOT/'docs/verification/paired-geometry-ancf-physical-pde-results.md';assert not report.exists()
    maximum=lambda key:max(max(a['metrics'][key]) for a in r['records'])
    rows=[]
    for a in r['records']:
        rows.append(f"| {a['mode']} | {a['amplitude']:g} | {a['n']} | {a['h']:g} | {max(a['metrics']['error']):.6g} | {max(a['metrics']['projection_field_energy_difference']):.4g} | {max(a['metrics']['relative_work_defect']):.4g} |")
    rational=[a for a in r['records'] if a['mode']=='rational' and a['n']==4];ratios=[a['ratio'] for a in r['accuracy'] if a['mode']=='rational' and a['n']==4 and a['field']=='error']
    text=f'''# R4.1：PG-ANCFの物理荷重平衡・反力・仕事

2026-11-04。**解析放物柱面・均質等方・線形静的の物理荷重経路を採用する。R4.1を終了する。** PG基底を変えず、独立な連続体の体積力と面力から平衡を解き、R3の投影診断との同値を確認した。一般曲面・材料／積層等のR4全体、主シェル最終採用、本体実装開始は未完了。

## 成立したこと

[事前仕様](paired-geometry-ancf-physical-pde-full-protocol.md)の正則／歪み、n2/n4、h=.01/.0001、標準／独立曲率の**16組・48右辺**を通過。各求解中は製造ひずみ／変位の評価器を例外に置き換え、物理外力から粗い初期場と全自由度補正を独立に求めた。保存済みの投影解は求解後の比較にだけ使用した。[全結果](paired-geometry-ancf-physical-pde-batch-v3-result.json)。

参照は既存のd₀=(−κx,0,1)の柱面で、実際の法線厚さはh√(1+κ²x²)。hは参照写像の厚さパラメータである。一定の実法線厚さを指定する幾何モデルは、次の一般参照面の検査へ残る。

| 検査量 | 最大値 | 閾値 |
|---|---:|---:|
| 保存投影解との場の差：参照エネルギーで正規化 | {maximum('projection_field_energy_difference'):.6g} | 1e−8 |
| 保存投影解との位置場差：参照L²位置で正規化 | {maximum('projection_field_position_difference'):.6g} | 1e−8 |
| 元点の自由平衡残差 | {maximum('stationarity'):.6g} | 1e−8 |
| 外力仕事と弾性仕事の相対差 | {maximum('relative_work_defect'):.6g} | 1e−8 |
| 合力・モーメントの釣合い | {maximum('relative_force_moment_balance'):.6g} | 1e−8 |
| 計算反力と連続体反力の合力・モーメント差 | {maximum('relative_reaction_resultant_difference'):.6g} | 1e−8 |
| 別57点の固定端変位 | {max(a['metrics']['fixed_boundary'] for a in r['records']):.6g} | 1e−9 |

独立曲率の最細n4エネルギー相対誤差は最大{max(max(a['metrics']['error']) for a in rational):.9g}、薄肉増幅は最大{max(ratios):.9g}倍。1%／2倍を維持。数値床以下の標準場の比から新しい精度・収束率を主張しない。格子減少、厚さ級数の独立差、資源対応前後の同じ物理場は[最終監査](paired-geometry-ancf-physical-pde-audit.json)で別に検査する。

| 荷重群 | 歪み振幅 | n | h | 最大エネルギー誤差 | 投影解との場の差 | 仕事差 |
|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

JSONのelastic_energyは**弾性仕事a(u,u)=2U**の値。外力仕事と比較する量であり、ポテンシャルUそのものではない。支持反力の一般化成分の一致を全成分で要求したわけではなく、物理合力・モーメントと固定端仕事の整合を検査した。

## 荷重組立と処理境界

右辺はb=−divσ、t=σnから作り、Kと既知変位の積を外力として与えていない。体積力と上下面力のO(h)項が相殺してO(h³)の曲げ荷重となるため、厚さモーメントの係数を荷重側で先にまとめた。側面は別に全て積分し、固定側のσnは反力参照として分離した。

[80桁の独立な直接積分](paired-geometry-ancf-load-moment-audit.json)と閉じたモーメントの差は最大5.18e−16。求積・粗い仕事・反力回復、測度・単位・適用面・失敗保証は[物理荷重の内部契約](../paired-geometry-ancf-physical-load-contract.md)にまとめた。汎用の任意荷重積分器や追従荷重の非線形受入まで成立したという判断ではない。

## 資源停止と同じ作用素の効率化

[元系列](paired-geometry-ancf-physical-pde-result.json)は12組を通過後、歪みn4・厚い側の独立曲率で数値条件を全通過したが、311.625秒で1配置300秒の枠を超えた。その不合格を保持する。[2048点経路の準備後の中断](paired-geometry-ancf-physical-pde-stream-interruption.json)も記録した。

[拘束後QR](paired-geometry-ancf-physical-pde-constrained-protocol.md)では、同じ支持写像Zを先に適用しB(S)Z=B(SZ)から全747許容成分を組み立て、荷重で不要な導関数の基底変換を省いた。作用素／位置は一致したが318.579秒で時間未達となり、[その結果](paired-geometry-ancf-physical-pde-constrained-result.json)も保持した。

[参照キャッシュ経路](paired-geometry-ancf-physical-pde-cached-protocol.md)では、基底の配列ハッシュと参照・点をキーとした読取専用キャッシュで、同じ形状関数値と粗い基底を共有した。荷重や状態に依存する値は共有しない。点・重み・材料・荷重・許容未知量・閾値を変えず、作用素／位置／基底の元評価との別点一致を1e−12で要求し、元点の全許容自由度平衡も再確認した。参照行列の保存再利用時も元の準備時間を資源判定へ加算する。

キャッシュだけの経路も312.157秒で時間未達となり、[その結果](paired-geometry-ancf-physical-pde-cached-result.json)を保持。[中面評価経路](paired-geometry-ancf-physical-pde-plane-protocol.md)では同じ面内点の基底・導関数を一度評価して厚さ多項式で展開し、N・物理勾配・Bを元評価と別点で照合した。全積分点の元Jacobianと同じ全許容未知量を保持する。

中面評価の経路でも2荷重群を別々に観測した配置全体は376.235秒となり、[時間未達](paired-geometry-ancf-physical-pde-plane-result.json)を保持した。[6右辺の一括観測](paired-geometry-ancf-physical-pde-batch-protocol.md)では、異なる物理右辺を別列のまま同時に解き、同じN/Bを一度評価した。3列群へ戻して全判定を保存している。接続時の名前解決漏れによる[v1停止](paired-geometry-ancf-physical-pde-batch-result.json)・[v2停止](paired-geometry-ancf-physical-pde-batch-v2-result.json)も残し、修正したv3を判定対象とする。

通過済み12組を再利用し、残4組を同じ基準で再計算。残系列の今回の実行{r['seconds']:.1f}秒、ピーク{r['peak_bytes']/1024**3:.3f}GiB。保存済み行列の元の準備時間を加えた各残幾何配置＋6右辺は最大{max(a['configuration_elapsed'] for a in r['records'] if a.get('execution_route')=='six_rhs_shared_geometry'):.1f}秒で、300秒・2GiB枠内。これは試作の資源検査であり製品性能の受入ではない。

## 次に閉じる条件

次はR4の一般参照面：参照Ritz持上げ、必要二階微分、実際の法線厚さ、チャートと全トレースを一組で定義・照合する。柱面のκを点ごとに差し替えるだけの式で一般面を採用しない。その後、既に採用した材料点契約を要素へ接続し、積層／初期ひずみ、非線形・慣性／固有特性・再開、代表規模の容量を同じSへ対応させる。R4では採用を左右する代表照合と実装契約を閉じ、製品の網羅検証は実装と対応させる。

[全体ロードマップ](../roadmap.md)、[判断と出典](paired-geometry-ancf-physical-pde-disposition.json)。主シェルの最終採用はR5で判断する。
'''
    report.write_text(text,encoding='utf-8')
    message='PG-ANCF候補のR4.1を終了。独立な体積力・面力による16組48右辺の平衡、保存投影解との場の一致、反力・仕事・合力／モーメントが成立。物理荷重の内部I/Fと失敗保証を定義。次は一般参照面のRitz持上げ・二階微分・実法線厚さ／トレース、その後材料・慣性等との整合。R4全体・主シェル最終採用・本体実装開始は未完了。'
    for name,prefix,rel in [('README.md','2026-11-04：','docs/'),('docs/preimplementation-plan.md','現在の進行（2026-11-04）：',''),('docs/structural-model-adoption-review.md','現在の到達点（2026-11-04）：',''),('docs/requirements-evidence-map.md','E-SHELLの現在の追加証拠（2026-11-04）：','')]:
        p=ROOT/name;lines=p.read_text(encoding='utf-8-sig').splitlines();hits=[i for i,s in enumerate(lines) if s.startswith(prefix)];assert len(hits)==1;lines[hits[0]]=prefix+message+f'[R4.1判定]({rel}verification/paired-geometry-ancf-physical-pde-results.md)、[現在地]({rel}roadmap.md)。';p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'docs/design-decisions.md';lines=p.read_text(encoding='utf-8-sig').splitlines()
    for i,line in enumerate(lines):
        if line.startswith('| ST-01 |') or line.startswith('| ST-04 |'):
            cells=line.split('|');cells[-2]=' P2・'+message+'（[R4.1](verification/paired-geometry-ancf-physical-pde-results.md)、[荷重契約](paired-geometry-ancf-physical-load-contract.md)、[計画](roadmap.md)） ';lines[i]='|'.join(cells)
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'docs/roadmap.md';lines=p.read_text(encoding='utf-8-sig').splitlines()
    for i,line in enumerate(lines):
        if line.startswith('**現在地は'):lines[i]='**現在地は1.1主シェル、PG-ANCF候補のR4を進行中。R4.1の柱面・均質等方・線形静的の物理荷重経路を閉じた。** [16組48右辺・反力／仕事の判定](verification/paired-geometry-ancf-physical-pde-results.md)、[物理荷重I/F](paired-geometry-ancf-physical-load-contract.md)。次は一般参照面のRitz持上げ・二階微分・実法線厚さ／全トレース、その後材料・慣性・費用を同じSへ対応させ、R5で最終採用を判断する。R2/R3候補成立と旧方式の不採用は維持。主シェル最終採用・本体実装開始は未承認。'
        elif line.startswith('| **R4 実務'):lines[i]='| **R4 実務に必要な力学との整合：進行中** | **R4.1の柱面物理荷重・平衡・反力・仕事を終了**。[判定](verification/paired-geometry-ancf-physical-pde-results.md)。次は一般面のRitz持上げ・二階微分・実法線厚さ／トレース。その後、材料／積層・初期ひずみ、非線形・慣性／固有特性・再開、代表規模の容量を同じSへ接続 | 代表照合、適用／拒否条件、主要I/Fと不変条件、実装後へ残す検証を対応させる。R4全体は未完了。一般面や費用が未達なら最終採用しない |'
        elif line.startswith('| **1 構造モデル採用**'):lines[i]='| **1 構造モデル採用** | **進行中・主シェル未採用**。PGの候補選別に加えR4.1の物理荷重経路を終了 | 一般面の参照Ritz・厚さ／トレース→材料・慣性・費用の整合→R5採用。その後、梁・接続等を統合 | ST-01〜10/CV-03の採用式・適用条件・必要微分・接続・独立参照と検証計画が揃う |'
        elif line.startswith('| 1.1 主シェル'):lines[i]='| 1.1 主シェルの補間・ロッキング・求積・幾何 | **PG候補のR2/R3とR4.1物理荷重を終了**。[反力・仕事と独立平衡](verification/paired-geometry-ancf-physical-pde-results.md)が成立 | 一般面・実法線厚さ／トレース、材料・非線形・慣性・容量とR5最終採用は未完了 |'
        elif line.startswith('R1〜R5は'):lines[i]='R1〜R5は実装前6作業の「1」の中で主シェルを閉じる作業。PG候補のR2/R3を終え、R4.1で物理荷重による平衡と投影診断の同値も確認した。次は一般参照面を扱う。R4全体とR5は未完了で、項目数を完了率や採用日に換算しない。'
        elif line.startswith('| **R3 決定的'):lines[i]=line.replace('物理PDE・一般面・材料等の最終採用は未完了','物理PDEはR4.1で別に検査。一般面・材料等の最終採用は未完了')
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    paths=[path,report,ROOT/'docs/paired-geometry-ancf-physical-load-contract.md',ROOT/'docs/verification/paired-geometry-ancf-load-moment-audit.json',Path(__file__),ROOT/'docs/verification/paired-geometry-ancf-candidate-audit.json']
    disposition=ROOT/'docs/verification/paired-geometry-ancf-physical-pde-disposition.json';assert not disposition.exists();disposition.write_text(json.dumps(dict(R4_1_bounded_physical_load_path_accepted=True,R4_complete=False,main_shell_adopted=False,product_implementation_started=False,scope='parabolic cylinder; homogeneous isotropic linear statics; body and all applicable surface forces',energy_quantity='a(u,u)=2U',next='general reference surface, Ritz lift, second derivatives, physical normal thickness and complete traces',sources={p.relative_to(ROOT).as_posix():sha(p) for p in paths}),indent=2),encoding='utf-8');print('R4.1 report, decision and roadmap updated')

if __name__=='__main__':main()
