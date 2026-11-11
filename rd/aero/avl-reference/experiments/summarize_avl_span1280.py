"""Publish result tables only from completed, audited numerical artifacts."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import json
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
def main():
    raw=json.loads((F/'avl-span1280-fine-ks-runs.json').read_text());obs=json.loads((F/'avl-span1280-observation-result.json').read_text());audit=json.loads((F/'avl-span1280-fine-artifact-audit.json').read_text())
    assert raw['completed'] and obs['methods_criteria_pass'] and audit['all_criteria_pass']
    def percent(v):return '未評価' if v is None else f'{100*v:.6f}%'
    names={'tau0':'τ帯平均','tau_p1':'τ帯P1','k0':'K帯平均','j_line':'有限辺モーメント密度','j_point':'元点モーメント密度'}
    lines=['# N1280・同一基準点の保存観測：実施結果','',
        '2026-10-31。[細分化事前仕様](avl-span1280-fine-experiment.md) に従い、対称静的2翼の物理全翼N1280までを計算した。参照全体の採用は保留する。','',
        '## 1. 初回失敗と別版の修正','',
        '初回はFW-V01・M3のLOAD段階でamake.fの局所KSMAX=500が不足して停止した。[失敗記録](avl-span1280-fine-runs.json) と元版は保持した。Fortran STOPの終了コードは0だったが、console検査で拒否し、解なしの結果を後処理へ渡していない。前回のNVX／NVMAX／NSMAXの調査は、この局所容量を見落としていた。','',
        '[追加事前仕様](avl-span1280-ksmax-addendum.md) により別コピーでKSMAXだけを1000へ拡張し、amake.oを再コンパイル、他40オブジェクトを親とハッシュ一致で再利用した。[8入力の読込み](avl-span1280-load-result.json)、既存8ケース・66出力のバイト一致と [修正版監査](avl-span1280-capacity-ks-artifact-audit.json) が成立してから再実行した。格子・核・求解法・判定しきい値は変更していない。','',
        '通常consoleやFEのI4整数欄では15360を表示できない。[表示幅の追加仕様](avl-span1280-load-diagnostics-addendum.md) に従って制限を記録し、未知数・帯数・力と場の評価はUH metaおよび高精度行／力ファイルで確認した。通常出力の書式変更を今回の修正に混ぜていない。修正版容量回帰のマニフェストに残った親ビルド資源ハッシュはそのまま保持し、実際の修正版ビルド資源は修正版監査で別に検証・ハッシュ化した。','',
        '## 2. 計算と資源','',
        '[修正版準備](avl-span1280-fine-ks-preparation.json)、[8計算の実行記録](avl-span1280-fine-ks-runs.json)、[停止監視の試験](avl-span1280-watchdog-result.json) を保持する。通常／拒否の監視分岐と実際の子プロセス停止を検証後に実行した。','',
        '| 最細ケース | 実未知数 | 物理全翼パネル | 時間 | ピークRSS bytes |','|---|---:|---:|---:|---:|']
    for r in raw['runs']:
        if r['m']==24:lines.append(f'| {r["case"]}・M24N1280 | {r["native_unknowns"]} | {2*r["native_unknowns"]} | {r["seconds"]:.2f}秒 | {r["peak_resident_bytes"]:,} |')
    lines+=['',f'8件の合計は{raw["seconds"]:.2f}秒、元出力の保存量は{raw["total_bytes"]:,} bytes。各1200秒・計2400秒・群RSS16GiB・開始有効空き20GiB・残量2GiB・swap増加256MiB・出力512MiBの条件内だった。仮想割当量と常駐メモリ、major faultとswapは別記録とした。', '',
        '## 3. 保存・残差・場間差上限','',
        f'[観測結果](avl-span1280-observation-result.json) は同じKSMAX修正版のN160／320／640回帰とN1280から14場を構成した。N1280の半翼640帯・最細15360未知数をUH metaと行数で確認し、物理全翼へ別ファイルとして展開した。元残差は最大相対{max(r["residual_relative"] for r in obs["meshes"]):.6e}。処理検査{len(obs["checks"])}条件は全て成立した。','',
        '計測前後の出力不変、全域合力／元点モーメント、τ/K恒等式、固定帯の独立集計、帯1次モーメントを確認した。12組の二つの保存場間の上限も成立した。これらは未知の連続解への誤差上界ではなく、非剛体のANCF荷重伝達の証拠でもない。','',
        f'後処理は{obs["seconds"]:.2f}秒、ピーク{obs["peak_process_working_bytes"]:,} bytesで120秒・1GiB以内。[1翼目の先行評価](avl-span1280-first-wing-result.json) は [別チェックポイント](avl-span1280-first-wing-checkpoint.json) に結び付けて保持し、最終の同じ翼の系列・予算と一致することも監査した。','',
        '## 4. 同じ判定条件による4系列','',
        '両翼それぞれでM24・N160／320／640／1280と、N1280・M3／6／12／24を評価した。差の減少・最終次数0.5〜6・最終2次数の変化0.5以下・両方向余弦0.9以上・丸め床は従来どおり。','',
        '| 翼・系列 | チャンネル | 状態 | 最終2次数 | 方向余弦2値 | 最終相対差 |','|---|---|---|---|---|---:|']
    for s in obs['series']:
        for key,c in s['channels'].items():
            orders='—' if c['orders'] is None else ' → '.join(f'{v:.5f}' for v in c['orders'])
            cosines=' / '.join('—' if v is None else f'{v:.5f}' for v in c['difference_cosines'])
            status={'asymptotic_candidate':'漸近候補','not_asymptotic':'未達','roundoff_candidate':'丸め床候補'}.get(c['status'],c['status'])
            lines.append(f'| {s["case"]}・{s["label"]} | {names[key]} | {status} | {orders} | {cosines} | {percent(c["relative_last_difference"])} |')
    lines+=['','## 5. M24N1280での仮定付き予算','',
        '| 翼 | 揚力 | 帯分布 | 元点モーメント | 全目標 |','|---|---:|---:|---:|---|']
    for b in obs['budgets']:lines.append(f'| {b["case"]} | {percent(b["conditional_lift_relative"])} | {percent(b["conditional_distribution_relative"])} | {percent(b["conditional_moment_relative"])} | {"条件付きで予算内" if b["conditional_targets_met"] else "未達／未評価あり"} |')
    lines+=['',
        '予算は同じ基準点の弦／スパン推定を線形和する定義で、端点移送を使わない。今回はスパン推定が成立せず、合算値は未評価とした。目標は揚力・元点モーメント0.3%、帯分布0.5%のまま。必要な方向が候補でなければ対応値はnullのまま保持する。','',
        '外挿は漸近尾部の仮定に依存し、方向間の混合尾部は保証していない。固定20帯の重み付き観測の推定を、任意点の相対誤差や未知の連続解への保証に読み替えない。','',
        '## 6. 判定と次の作業','',
        '[成果物監査](avl-span1280-fine-artifact-audit.json) は、初回失敗、修正版ゲート、8入力・実行体・ソース・出力・派生場・資源・条件付き予算の追跡を検証した。参照全体の採用、UVLAS対参照の一致、実験妥当性、非対称／動的連成、J1および製品数値コア開始の判定は更新しない。','',
        '次は、この表で未達の観測量と誤差仮定を整理し、既存の保存観測差の符号・方向・相殺から次数変動を診断する。参照として限定使用できる出力量・格子・不確かさの扱いは、その証拠に基づき明記する。N640で条件付き推定が得られていた量も、新しい窓で必要な条件を失った場合は、今回の基準点で未評価として扱う。','',
        '同じ経路の次のN2560・M24は半翼30720未知数・1280帯となり、現容量と局所容量を超える。NVX32000で同じ7枚相当を確保すれば57,344,000,000 bytes（約53.41GiB）で、今回のWSL総メモリ約23.38GiBを超える。実行全体の上界ではなく配列だけの値であり、そのままの追加細分化は実行できない。細分化を選ぶ場合には、メモリ／求解経路・時間予算・回帰と独立した判定根拠を先に定義する。','']
    (F/'avl-span1280-results.md').write_text('\n'.join(lines),encoding='utf-8')
    print('Wrote audited result tables')
if __name__=='__main__':main()
