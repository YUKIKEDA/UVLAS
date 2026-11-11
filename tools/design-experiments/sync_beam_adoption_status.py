"""Synchronize six mutable project entry points after immutable beam audit."""
from pathlib import Path
import hashlib,json,re
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=ROOT/'docs/verification/beam-adoption-audit.json';a=json.loads(ap.read_text());assert a['passed'] and a['beam_1_3_complete']
    for label in ['sources','artifacts']:
        for name,h in a[label].items():assert sha(ROOT/name)==h,name
    common='1.3を完了し、SE(3)幾何学的厳密せん断変形梁を採用。軸・せん断・ねじり・曲げ、12接自由度、連成断面、初期ひずみ、構成依存質量／慣性、シェルへの点・フレーム・J/Hと転置力を固定した。20群と25入力負対照が成立。必要な状態精度と拘束反り等の適用境界を明示。主シェル1.1、材料1.2、梁1.3は採用済みで、現在地は1.4接続・支持・質量・偏心慣性。構造全体・本体実装開始は未完了。'
    fronts=['README.md','docs/preimplementation-plan.md','docs/structural-model-adoption-review.md','docs/requirements-evidence-map.md']
    for name in fronts:
        p=ROOT/name;lines=p.read_text(encoding='utf-8').splitlines();count=0
        for i,line in enumerate(lines):
            marker='1.2材料・積層・等価入力の数値仕様と主シェル結合を採用。'
            if marker in line:
                prefix=line.split(marker)[0];base='docs/' if name=='README.md' else ''
                lines[i]=prefix+common+f' [梁採用判定]({base}verification/beam-adoption-results.md)、[梁採用仕様]({base}beam-adopted-specification.md)、[現在地]({base}roadmap.md)。';count+=1
            if name=='README.md' and line.startswith('**現在地：'):lines[i]='**現在地：実装前6作業の「1. 構造モデル採用」→1.4 接続・支持・質量・偏心慣性。** 主シェル1.1、材料1.2、梁1.3は採用済み。構造審査全体は未完了。C++23の両OSビルド・CLI試験の基盤はあるが、製品の構造・空力・連成計算本体は未実装。'
        assert count==1,(name,count);p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'docs/design-decisions.md';lines=p.read_text(encoding='utf-8').splitlines()
    states={
        'ST-06':'P2・1.3完了。BM-D01〜09として2節点SE(3)せん断変形梁、位置／回転と12接自由度、全連成断面、一般化初期ひずみ、構成依存質量／慣性、点・フレーム・J/H・転置力を採用。20群・25入力負対照。拘束反りは別経路、必要精度と投影適格性を4へ渡す',
        'ST-05':'P2・PG-ANCFの一定整合質量とSE(3)梁の構成依存M(q)/Christoffelを採用。偏心断面の物理慣性、固有対、運動エネルギーを独立照合。集中質量・一般偏心剛体・支持／接続との統合は1.4',
        'ST-08':'P2・主シェルの全トレース、梁の点／フレームとJ/H・転置力、既存厚さ座標変換を入力として固定。現在1.4で実際の梁―シェル、折れ・交差・剛体・指定運動、独立拘束・反力・冗長／不整合を統合する。接合全体の採用は未完了',
        'ST-10':'P2/P3・採用シェル／材料とSE(3)梁の同じエネルギーからの内力・材料＋幾何接線・回復境界を採用。独立高精度微分と物理参照が成立。梁の動的残差の高階微分は3／4、接合／断面力の全体回復は1.4／1.5で統合'}
    for i,line in enumerate(lines):
        for key,text in states.items():
            if line.startswith('| '+key+' |'):
                cells=line.split('|');cells[-2]=' '+text+'（[梁仕様](beam-adopted-specification.md)、[1.3判定](verification/beam-adoption-results.md)） ';lines[i]='|'.join(cells)
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    # The integration table is a live status table, distinct from frozen tests.
    p=ROOT/'docs/structural-model-adoption-review.md';lines=p.read_text(encoding='utf-8').splitlines()
    for i,line in enumerate(lines):
        if line.startswith('| ST-01/04 '):lines[i]='| ST-01/04 主シェル・求積・安定性 | [PG-ANCF採用仕様](main-shell-adopted-specification.md)、[R4整合](verification/paired-reference-r4-results.md)、[R5採用](verification/paired-shell-r5-results.md) | **1.1の採用判断は完了**。製品全実務受入・具体数値依存は後続 | 採用式、必要微分、適用／拒否、全トレース、求積と保存証拠を固定済み |'
        elif line.startswith('| ST-02 '):lines[i]='| ST-02 参照配置・幾何 | [参照係数空間](paired-reference-coefficient-definition.md)、[採用仕様](main-shell-adopted-specification.md)、非二次・後退／テーパの証拠 | 多部材の参照面対応、外部格子変換、幾何誤差と実務品質条件は1.4／4／5 | 入力から参照写像への対応、品質診断と負対照、変換誤差予算 |'
        elif line.startswith('| ST-03/CV-03 '):lines[i]='| ST-03/CV-03 材料・積層・初期ひずみ | [材料採用仕様](material-shell-adopted-specification.md)、[1.2検証・採否](verification/material-shell-adoption-results.md) | **1.2の数値仕様・主シェル結合は完了**。梁は別の一般化ひずみ／断面則へ縮約。実材料の用途別受入は5 | 全3D入力・等価則、層／軸／初期ひずみ、全接線・質量・両側回復の境界を固定済み |'
        elif line.startswith('| ST-05/07 '):lines[i]='| ST-05/07 質量・偏心慣性 | シェル一定Mと梁M(q)、偏心断面、[物理空間慣性の検査](verification/beam-se3-core-v2-result.json) | **1.4の対象**。集中質量・偏心剛体・接続の全運動学と慣性を統合 | 運動エネルギー、一般化慣性、剛体特性・転置仕事・必要微分 |'
        elif line.startswith('| ST-06 '):lines[i]='| ST-06 梁 | [SE(3)梁採用仕様](beam-adopted-specification.md)、[1.3検証・採否](verification/beam-adoption-results.md) | **1.3の採用判断は完了**。実接合は1.4、動的方式／具体精度は3／4、実務範囲は5 | 軸・せん断・ねじり・曲げ、断面／反り境界、DOF・接線・質量、シェル受渡しを固定済み |'
        elif line.startswith('| ST-08 '):lines[i]='| ST-08 界面・支持・剛体結合 | 全トレース、梁の物理点／フレームのJ/Hと転置力、[厚さ変換](thickness-coordinate-interface.md) | **現在地1.4**。梁―シェル、折れ・交差、非一致、剛体、指定運動、独立拘束・反力 | C/J/H・時間微分、許容仮想変位、冗長・不整合の拒否と組合せ物理問題 |'
        elif line.startswith('| ST-10 '):lines[i]='| ST-10 内力・接線・結果回復 | シェル・材料・梁の採用エネルギー、全接線、点／層／断面回復の境界と独立照合 | 接合・支持・荷重・減衰と全体回復を1.4／1.5、生成方式／具体型を4で統合 | 同じ状態の処理仕様・必要微分、力／モーメント保存と回復可否 |'
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'docs/roadmap.md';lines=p.read_text(encoding='utf-8').splitlines()
    for i,line in enumerate(lines):
        if line.startswith('**現在地は1.3梁の採用。'):lines[i]='**現在地は1.4接続・支持・質量・偏心慣性。主シェル1.1、材料1.2、梁1.3は採用済み。** '+common+' [梁採用判定](verification/beam-adoption-results.md)、[梁仕様と1.4への引継ぎ](beam-adopted-specification.md)。'
        elif line.startswith('R1〜R5は実装前6作業'):lines[i]='R1〜R5で主シェル1.1を閉じ、材料1.2と梁1.3も採用済み。次は1.4の実接続・支持・質量／偏心慣性を統合する。候補探索を根拠なく再開せず、項目数を製品完成率や採用日に換算しない。'
        elif line.startswith('- R4/R5と材料1.2の成立範囲'):lines[i]='- 主シェルR4/R5、材料1.2、梁1.3の証拠・仕様を再利用する。次の有限な比較は1.4の実接合・独立拘束・反力・慣性と組合せ物理問題に限定する。'
        elif line.startswith('| **1 構造モデル採用**'):lines[i]='| **1 構造モデル採用** | **進行中・主シェル1.1、材料1.2、梁1.3を採用済み** | **1.4接続・支持・質量／偏心慣性→1.5荷重・回復・構造統合** | ST-01〜10/CV-03の採用式・適用条件・必要微分・接続・独立参照と検証計画が揃う |'
        elif line.startswith('| **1.3 梁・シェル受渡し'):lines[i]='| 1.3 梁・シェル受渡し | **採用完了**。[SE(3)梁仕様](beam-adopted-specification.md)、20群・25入力負対照、独立解析／連続体参照が成立 | 実接合は1.4、時間・高階慣性微分は3、具体精度／投影と数値依存は4、実務受入は5へ条件を固定して引継ぎ |'
        elif line.startswith('| 1.4 接続・支持・質量'):lines[i]='| **1.4 接続・支持・質量・偏心慣性：現在地** | シェル全トレースと梁の点／フレーム・J/H・転置力、個別慣性の証拠は揃った | 梁―シェル・折れ・交差・厚さ段差・剛体・指定運動、独立拘束・反力・集中／偏心慣性を同じ運動学から閉じ、組合せ物理問題を比較 |'
        elif line.startswith('1.1と1.2の採用判断を閉じたため'):lines[i]='1.1〜1.3の採用判断を閉じたため、同じ主シェル・材料・梁の運動学から実際の接合網と慣性を統合できる。1.4では受渡しの仮想仕事だけで終わらず、参照適合・独立拘束・組合せ平衡／反力／固有特性まで確認する。'
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    for label in ['sources','artifacts']:
        for name,h in a[label].items():assert sha(ROOT/name)==h,name
    allfronts=fronts+['docs/design-decisions.md','docs/roadmap.md'];links=0
    for name in allfronts:
        p=ROOT/name;s=p.read_text(encoding='utf-8');assert '1.4' in s and 'beam-adopt' in s,name
        for target in re.findall(r'\]\(([^)]+)\)',s):
            if 'beam-adopt' in target or 'beam-se3' in target:
                assert (p.parent/target.split('#')[0]).exists(),(name,target);links+=1
    frozen=['docs/beam-adopted-specification.md','docs/verification/beam-adoption-disposition.json','docs/verification/beam-adoption-results.md','docs/verification/beam-adoption-audit.json']
    closure=dict(passed=True,beam_1_3_complete=True,material_1_2_complete=True,main_shell_adopted=True,structural_review_complete=False,product_implementation_start=False,
        current_work='1.4 connections supports mass and eccentric inertia',immutable={name:sha(ROOT/name) for name in frozen},status_sync_source={Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))},mutable_status_entries=allfronts,checked_status_links=links)
    p=ROOT/'docs/verification/beam-adoption-closure.json';assert not p.exists();p.write_text(json.dumps(closure,indent=2),encoding='utf-8');print(json.dumps(closure,indent=2))
if __name__=='__main__':main()
