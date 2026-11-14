"""Summarize independent curved consistency failures and fixed design constraints."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json
ROOT=rd_paths.ROOT;F=rd_paths.RECORDS
def read(n):return json.loads((F/n).read_text())
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    data=read('q9-curved-strain-screen-result.json');audit=read('q9-curved-strain-screen-audit.json');boundary=read('q9-curved-strain-boundary-result.json');ba=read('q9-curved-strain-boundary-audit.json')
    assert data['diagnostic_checks_pass'] and boundary['diagnostic_checks_pass'] and audit['all_criteria_pass'] and ba['all_criteria_pass']
    assert audit['parent_sha256']==digest(F/'q9-curved-strain-screen-result.json') and ba['parent_sha256']==digest(F/'q9-curved-strain-boundary-result.json')
    lines=['# 一般曲面の一定ひずみと節点仕事の審査結果','','2026-11-02。**T-S36-v1の素朴な3D共変拡張G-COVと、中心基底へ交換したC-COVは、一般曲面用には採用しない。** 平坦のT-S36-v1の限定比較通過は維持し、不合格の曲面拡張と区別する。','','## 定義した範囲','','[事前仕様](q9-curved-strain-screen-experiment.md) に従い、RECT/SKEW/TRAP/SADDLE/NORMAL-C、厚さ0.1／0.01、参照配置の回転、厚さ座標の再表現を組み合わせた40配置を、2作用素で検査した。荷重平衡や製品数値コアは実装していない。','','G-COVは実際の参照基底で全ひずみを共変に変換し、13／23だけを36点から補間して全テンソルを戻す。C-COVは中心の基底を全点で使う比較案。後者はMITC9iの全構成を実装したものではなく、その論文・要素一般を不適格とする判定ではない。','','## 一定場の値と一定応力の節点仕事','','下表は回転・再表現前の配置。値差は6一定物理ひずみの最大相対テンソル差、仕事差は6一定物理応力に対する全81仮想変位の節点力作用素の最大相対差。変位応答の誤差率ではない。いずれも事前基準1e−10で判定した。','','| 形状 | 厚さ | 案 | 一定場の値差 | 一定応力の仕事差 |','|---|---:|---|---:|---:|']
    for c in data['cases']:
        if c['rotated'] or c['reparameterized']:continue
        lines.append(f"| {c['geometry']} | {c['h']} | {c['scheme']} | {c['value_error']:.6e} | {c['work_error']:.6e} |")
    lines += ['', 'G-COVはSADDLEで一定場の値を保持しても、仕事差が約2.87%となった。NORMAL-Cでは値差約1.17%、仕事差約2.77%。C-COVは全形状で値を保持するが、SADDLEの仕事差は約0.875%、NORMAL-Cは約0.794%だった。値の再現を仕事の成立の代用にしない。','','全80記録の原始アフィン場、回転・厚さ再表現の共変性、同じ物理幾何の対応は事前条件内。7／9点の仕事差も全て事前基準内で、作用素の欠損を積分未収束へ付け替えない。成立した処理条件と不合格の採用ゲートを分けて保存した。','','[別の双対応力による監査](q9-curved-strain-screen-audit.json) は回転・再表現前の全10形状厚さ×2案で成立。共変成分を直接再構成する主経路と、保持応力・サンプルへ戻す双対寄与から節点力を作る別経路の差は最大約1.71e−15。元の一定応力節点力は、6境界面のσnを直接積分した結果とも約1.82e−14以内で一致した。人による独立審査ではない。','','## 辺上の欠損と共有トレース','','[追加境界仕様](q9-curved-strain-boundary-experiment.md) で4辺の一定場と、二要素・135仮想変位の接線－厚さトレースを検査した。下表のトレース差は作用素ノルム比であり、荷重応答の誤差率ではない。','','| 形状 | 厚さ | 案 | 辺上の一定場の値差 | 共有トレース差 |','|---|---:|---|---:|---:|']
    for c in boundary['cases']:
        if c['geometry'] not in ['SADDLE','NORMAL-C','TRAP']:continue
        lines.append(f"| {c['geometry']} | {c['h']} | {c['scheme']} | {c['boundary_constant_value_error']:.6e} | {c['shared_trace_error']:.6e} |")
    lines += ['', 'G-COVは共有トレースを保持するが、NORMAL-Cの辺上に約1.30%の一定場欠損を残す。C-COVは一定場の値を保持しても、台形・曲面の共有トレースで約0.24の相対作用素差となった。元の非還元トレースと参照接線は両側で一致し、入力の隙間・不適合を原因とする結果ではない。', '', '[別の6成分写像と点ごとの双対を用いた境界監査](q9-curved-strain-boundary-audit.json) でも全20組の判定が一致し、主経路との差は約2.24e−16以下だった。','','## 確定した判断','','- G-COVによるそのままの曲面拡張は不採用。一定物理ひずみと一定応力の節点仕事を同時に保持できない。','- C-COVへの中心基底の単純交換も不採用。値再現だけでは仕事・共有トレースを閉じられない。','- G-COVの辺上の値欠損、C-COVの共有辺欠損は、辺で0となる内部bubbleだけでは修復できない。','- [作用素の整合条件](../q9-curved-reduction-consistency-contract.md) を固定した。一定場の再現R P=P、一定応力の仕事P R=P、共有トレース、安定な無せん断空間を同時に課す。節点力は還元作用素の双対から作り、評価点の通常Bだけで代用しない。','- 次は幾何依存の辺還元・内部拡張とモーメント条件の可解性を先に定式化する。初期ひずみ・異方材を不合格作用素に追加して調整せず、荷重の大域計算も先行しない。','','旧P-S36の不採用、平坦T-S36の通過範囲、製品の要求と実装開始条件は維持する。主シェル最終採用・構造審査全体は未完了。','']
    path=F/'q9-curved-strain-screen-results.md';path.write_text('\n'.join(lines),encoding='utf-8')
    names=['q9-curved-strain-screen-result.json','q9-curved-strain-screen-audit.json','q9-curved-strain-boundary-result.json','q9-curved-strain-boundary-audit.json','q9-curved-strain-screen-results.md']
    result=dict(source_sha256=digest(__file__),evidence_sha256={n:digest(F/n) for n in names},screen_diagnostic_completed=True,independent_audits_pass=True,G_COV_accepted=False,C_COV_accepted=False,
        flat_T_S36_evidence_retained=True,next_work='geometry-dependent edge reduction and work-consistent extension feasibility',product_shell_adopted=False)
    (F/'q9-curved-strain-screen-disposition.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
