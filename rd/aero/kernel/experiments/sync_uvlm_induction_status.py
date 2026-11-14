"""Six live entrypoints: 2.1 complete, current 2.2, full UVLM still open."""
import sys,pathlib;sys.path.insert(0,str(next(p for p in pathlib.Path(__file__).resolve().parents if (p/'rd_paths.py').is_file())));import rd_paths
from pathlib import Path
import hashlib,json,re
ROOT=rd_paths.ROOT
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=ROOT/'rd/aero/kernel/records/uvlm-induction-adoption-audit.json';a=json.loads(ap.read_text());assert a['passed'] and a['uvlm_2_1_complete'] and not a['uvlm_review_complete'];dest=ROOT/'rd/aero/kernel/records/uvlm-induction-adoption-closure.json';assert not dest.exists()
    def verify():
        for group in ['sources','artifacts']:
            for name,h in a[group].items():assert sha(ROOT/name)==h,name
    verify();fronts=['README.md','rd/planning/preimplementation-plan.md','rd/planning/structural-model-adoption-review.md','rd/planning/requirements-evidence-map.md','docs/specs/design-decisions.md','rd/planning/preimplementation-history.md'];common='構造モデル採用1は完了。UVLMは2.1の有限3D W-P1誘導・近接／正則化自己評価・幾何微分を採用し、速度81列と幾何微分117列を個別照合、10群を検証。現在地は2.2の後流・組立・状態・荷重。2.3の有限翼／非定常の独立参照と予算、2.4のAE統合は未完了。UVLM全体・既定半径・自由後流更新・本体実装開始は未完了。'
    for name in fronts[:4]:
        p=ROOT/name;lines=p.read_text(encoding='utf-8').splitlines();count=0;base='docs/' if name=='README.md' else ''
        for i,line in enumerate(lines):
            if '1.5の荷重・減衰・支持・応力／断面力回復・時系列・状態規約を統合し' in line:
                lines[i]=line.split('1.5の荷重・減衰・支持・応力／断面力回復・時系列・状態規約を統合し')[0]+common+f' [2.1採用判断]({base}verification/uvlm-induction-adoption-results.md)、[誘導仕様]({base}uvlm-induction-adopted-specification.md)、[構造の終了判断]({base}verification/structural-adoption-results.md)、[現在地]({base}roadmap.md)。';count+=1
            if name=='README.md' and line.startswith('**現在地：'):lines[i]='**現在地：実装前6作業の「2. UVLM採用」→2.2 後流・組立・状態・荷重。** 構造1とUVLMの局所誘導2.1は完了。UVLM全体・連成／時間／安定性・具体数値基盤・実務受入・実装開始審査は未完了。'
        assert count==1,(name,count);p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'docs/specs/design-decisions.md';lines=p.read_text(encoding='utf-8').splitlines();states={
        'AE-01':'P2・dual-quarter-v1候補とW-P1-stream/P0-span-facet-v1の局所3D誘導を2.1として採用。幾何頂点／列×帯μを分離し、区分三角形の全基底を保持。全体の配置・安定ID・荷重との統合は現在2.2、最終採用は2.4',
        'AE-02':'P2・既存の限定非貫通／調和応答を保持。2.1で近接・自己・全基底の幾何微分を採用。一般3Dの非貫通組立とKutta／全荷重は2.2、独立有限翼精度は2.3、方式統合は2.4',
        'AE-03':'P2/P3・UI-D01〜06として同じ積分核の有限3D W-P1誘導・近接各基底・正則化自己・点／頂点／半径Jacobianを採用し2.1完了。10群、81速度列・117微分列。既定半径・成長・全域ID／交差・自由後流は2.2／2.3',
        'AE-04':'P2・上下／2軸ピッチのF/Mと4096分割の局所区間平均・空間変化の限定証拠を維持。高低保持と投影拒否を継承。有限3Dの誘導2.1は完了したが、全荷重／写像／幾何接線は現在2.2、独立精度は2.3',
        'AE-05':'P2・既存のW-P1放出・輸送・周期履歴・有限閉後流を入力に、正則化自己誘導と局所微分2.1を採用。現在2.2で核所有・物質ID・移流・初期化・反復内更新・確定／棄却を一組にする。自由後流更新の採用は未完了',
        'AE-06':'P2・局所共有辺・帯トレース不連続・共有頂点の鎖則を2.1で検証。複数揚力面／各後流の全相互組立、全域交差・近接と状態の適用診断は現在2.2。独立荷重参照は2.3'}
    for i,line in enumerate(lines):
        for key,value in states.items():
            if line.startswith('| '+key+' |'):
                cells=line.split('|');cells[-2]=' '+value+'（[誘導仕様](uvlm-induction-adopted-specification.md)、[2.1判定](verification/uvlm-induction-adoption-results.md)） ';lines[i]='|'.join(cells)
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    p=ROOT/'rd/planning/preimplementation-history.md';lines=p.read_text(encoding='utf-8').splitlines()
    for i,line in enumerate(lines):
        if line.startswith('**現在地は2 UVLM採用。'):lines[i]='**現在地は2 UVLM採用→2.2 後流・組立・状態・荷重。構造1と誘導2.1は完了。** '+common+' [2.1判定](verification/uvlm-induction-adoption-results.md)、[局所誘導仕様](uvlm-induction-adopted-specification.md)、[構造統合仕様](structural-adopted-specification.md)。'
        elif line.startswith('| **2 UVLM採用：現在地**'):lines[i]='| **2 UVLM採用：現在地2.2** | **2.1の有限3D誘導を採用済み**。既存の限定調和応答・局所平均圧力を保持。全体採用は未完了 | **2.2後流・組立・状態・荷重→2.3独立精度参照／予算→2.4 AE統合** | AE系の採用式・履歴・荷重と独立参照、数値誤差予算・拒否条件が揃う |'
        elif line.startswith('主シェル採用後は、構造1.2〜1.5を統合して1を閉じ'):lines[i]='構造1.1〜1.5の採用は終了。現在は2 UVLMを進め、その後3連成・時間・安定性、4数値基盤・主要I/F、5実務ケース・受入計画、6実装開始審査へ進む。既定の精度基準と適用範囲を維持する。'
    text='\n'.join(lines)+'\n';marker='## 現在地を構造審査の中で拡大する';assert marker in text
    block='''## UVLM採用の4ゲート

| ゲート | 状態・成立根拠 | 残作業と終了条件 |
|---|---|---|
| 2.1 有限3D誘導と局所微分 | **採用完了**。[10群・各速度81列／微分117列](verification/uvlm-induction-adoption-results.md)、[仕様](uvlm-induction-adopted-specification.md) | 指定核・半径の局所作用素として使用。既定半径・自由後流時間更新や有限翼精度の採用ではない |
| **2.2 後流・組立・状態・荷重：現在地** | 局所誘導、既存放出／履歴・荷重・幾何契約、採用PGシェルの点／フレームが入力 | 核所有・自己／共有ID・全域交差、全非貫通・放出／移流・初期化、試行／確定と有限幅の荷重／接線を同じ幾何で閉じる |
| 2.3 独立精度参照・誤差予算 | 幅一列の上下／ピッチ・固定区間圧力は限定成立。有限翼の受入参照不足は維持 | 定常／非定常、後退テーパ・有限幅・複数面と後流／コアの独立収束を審査。終了したAVL系列を既定で再開しない |
| 2.4 AE統合採用 | 未完了 | AE-01〜07の式・適用／拒否・必要微分・主要I/F・根拠を一組にし、UVLM採用2を閉じて3へ渡す |

4ゲートは審査の主対象を固定するための区切りであり、完了項目数を製品完成率へ換算しない。2.1だけで有限翼の荷重精度や自由後流の安定性を合格にしない。

'''
    text=text.replace(marker,block+marker);p.write_text(text,encoding='utf-8');verify();links=0
    for name in fronts:
        p=ROOT/name;s=p.read_text(encoding='utf-8');assert '2.2' in s and 'uvlm-induction' in s,name
        for target in re.findall(r'\]\(([^)]+)\)',s):
            if 'uvlm-induction' in target:assert (p.parent/target).exists(),(name,target);links+=1
    frozen=['docs/specs/aero/uvlm-induction-adopted-specification.md','rd/aero/kernel/records/uvlm-induction-adoption-disposition.json','rd/aero/kernel/records/uvlm-induction-adoption-results.md','rd/aero/kernel/records/uvlm-induction-adoption-audit.json'];closure=dict(passed=True,uvlm_2_1_complete=True,uvlm_review_complete=False,structural_review_complete=True,product_implementation_start=False,current_work='2.2 wake assembly state topology and loads',immutable={name:sha(ROOT/name) for name in frozen},status_sync_source={Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))},mutable_status_entries=fronts,checked_status_links=links);dest.write_text(json.dumps(closure,indent=2),encoding='utf-8');print(json.dumps(closure,indent=2))
if __name__=='__main__':main()
