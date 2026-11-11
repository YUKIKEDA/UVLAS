from pathlib import Path
import json,hashlib,re
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    audit_path=ROOT/'docs/verification/implementation-readiness-audit.json';audit=json.loads(audit_path.read_text());assert audit['passed'];assert audit['implementation_start_allowed_packages']==['IM-00'];parent=json.loads((ROOT/'docs/verification/numerical-foundation-adoption-closure.json').read_text());fronts=parent['mutable_status_entries']
    current='2026-11-06現在：実装前1〜4の採用設計、5の受入計画、6の最初の依存閉包審査を終了。IM-00の製品実装に着手可能。本体は未着手で、S10実験認定・製品規模性能・全実務受入は未完了。'
    def block(base):return '<!-- active-uvlm-status:start -->\n2026-11-06・最新：**実装前1〜4の採用設計、5の受入計画、6の最初の実装パッケージ審査を完了。次はIM-00：数値表現・求解・状態契約の製品実装。** 両OSの数値／Python接続を検証し、19ケース・37製品要求を対応付け、1072ソース・687成果物を監査。**S10実験参照の認定、実PG全規模・全自由後流の製品性能、全実務受入・リリースは未完了。本体コードは未着手。**\n\n[開始判定・実装順と残リスク]('+base+'implementation-start-review.md)、[受入計画]('+base+'practical-acceptance-plan.md)、[ケース台帳]('+base+'verification/practical-acceptance-case-catalog.json)。\n<!-- active-uvlm-status:end -->'
    for name in fronts:
        p=ROOT/name;s=p.read_text(encoding='utf-8');base='docs/' if name=='README.md' else ''
        s=re.sub(r'<!-- active-uvlm-status:start -->[\s\S]*?<!-- active-uvlm-status:end -->',lambda _:block(base),s,count=1)
        s=re.sub(r'^2026-11-05現在：構造1、UVLM2、J1共同決定3[^\n]*',current,s,flags=re.M)
        if name=='README.md':
            s='# UVLAS\n\n実務で利用できる品質のUVLM–ANCFシェル連成ソルバーを開発するプロジェクト。\n\n'+block('docs/')+'\n\n## 現在地\n\n方式・離散式・主要I/F・不変条件・参照・受入計画を整理し、最初の製品実装IM-00の開始条件を満たした。現在のsrcはビルド情報／CLI基盤のみ。設計実験のPython／C++コードは製品ソルバーではない。\n\n[全体ロードマップ](docs/roadmap.md)、[開始審査とIM-00〜08](docs/implementation-start-review.md)から作業と終了条件を確認できる。実験資料の不足、性能、精度消費、外部格子変換等の残リスクは開始審査R-01〜07で管理する。\n\n## 目標\n\n低速・非圧縮・付着流を対象に、複数揚力面と後流の相互干渉、等方材・積層複合材のANCFシェルと梁を組み合わせた構造、集中質量・偏心慣性・剛体結合、静的双方向連成、非線形過渡応答、大変形平衡周りのフラッタ予測を扱う。CLIとPython APIを提供し、WindowsとLinuxに対応する。\n\n代表性能目標は空力約5,000パネル・構造約10万自由度、単一PCでフラッタ速度探索全体12時間以内。主測定プロファイルは64GB搭載・プロセス48GiB。製品での達成実績はまだない。精度・適用範囲・測定条件は[受入計画](docs/practical-acceptance-plan.md)に定める。\n\n外部格子はNastran系の限定サブセットとGmsh。対応形式・カード・変換・拒否・再開は[入出力仕様](docs/io-execution-adoption-specification.md)に定める。途中確定状態と探索からの再開は必須。\n\n内部マイルストーンと公開リリースを区別する。初回の実務適用に必要な品質・機能を満たすまで部分公開しない。ソース、両OSのCLI、Python wheelを提供する計画。自身のソースはApache-2.0、第三者コード・データの条件は別に管理する。\n\n## 採用仕様と根拠\n\n- [製品要求](docs/product-requirements.md)、[精度・証拠条件](docs/acceptance-criteria.md)\n- [構造統合](docs/structural-adopted-specification.md)、[主シェルPG-ANCF](docs/main-shell-adopted-specification.md)\n- [UVLM統合](docs/uvlm-integrated-specification.md)、[連成・時間・安定性](docs/coupled-adopted-specification.md)\n- [初期・荷重微分の最終規約](docs/j1-initialization-output-derivative-addendum.md)\n- [数値基盤・型・I/F](docs/numerical-foundation-adoption-specification.md)、[入出力・実行・再開](docs/io-execution-adoption-specification.md)\n- [19ケースの受入計画](docs/practical-acceptance-plan.md)、[実装開始判定](docs/implementation-start-review.md)\n- [採用・不採用の経緯](docs/structural-model-adoption-review.md)、[設計台帳](docs/design-decisions.md)、[実装前計画](docs/preimplementation-plan.md)\n\n## ビルド基盤\n\nC++23。現在のビルド／CLI基盤の実行方法は[ビルド手順](docs/build.md)、両OSの確認は[C++23検証](docs/verification/cpp23-2026-10-29.md)。数値依存の限定試作は[4の終了判断](docs/verification/numerical-foundation-adoption-completion.md)。採用証拠の参照元は凍結し、製品へ移植した結果は別に検証する。\n'
        elif name=='docs/roadmap.md':
            a=s.index('## いま、どこにいるか');b=s.index('## 今後の計画',a)
            s=s[:a]+'## いま、どこにいるか\n\n**実装前準備の最初の依存閉包審査を終え、IM-00から製品実装に着手できる。** [開始判定・IM-00〜08・R-01〜07](implementation-start-review.md)。方式採用1〜4、受入計画5、最初の実装開始審査6は完了。本体は未着手。実務・実験・性能の合格は実装後の必須gateとして残る。\n\n'+s[b:]
            s=s.replace('2026-11-05更新。','2026-11-06更新。',1).replace('## 今後の計画（2026-11-05）','## 今後の計画（2026-11-06）')
            s=s.replace('F[5 実務ケース・受入計画：現在地]','F[5 受入計画：完了]').replace('G[6 実装開始審査]','G[6 最初の開始審査：完了]').replace('H[計算本体の実装・統合検証]','H[次：IM-00から製品実装]')
            s=s.replace('現在は5実務ケース・受入計画、その後6実装開始審査へ進む。','5実務ケース・受入計画と6の最初の開始審査も終え、次はIM-00から製品実装へ進む。')
            lines=s.splitlines()
            for i,line in enumerate(lines):
                if line.startswith('| **5 実務'):lines[i]='| 5 実務ケース・受入計画 | **計画採用完了**。A5-01〜07、19ケース・37要求、公開S10の取得／単位／不足を記録 | S10実験参照認定・実務検証結果・製品性能は実装後の必須gate | [5の終了判断](verification/practical-acceptance-plan-completion.md)。資料不足・未認定を合格にしない |'
                elif line.startswith('| 6 実装開始審査 |'):lines[i]='| 6 実装開始審査 | **最初の依存閉包は合格**。IM-00の開始条件8項目、IM-00〜08とR-01〜07を固定 | 次はIM-00の製品実装。後続は依存gateで進める | [開始判定](implementation-start-review.md)、[監査](verification/implementation-readiness-audit.json)。本体は未着手 |'
            s='\n'.join(lines)+'\n'
        elif name=='docs/preimplementation-plan.md':
            s=s.replace('版：1.2、2026-11-04。現状判定：計算本体の実装開始条件は未達。','版：1.3、2026-11-06。現状判定：最初の製品実装IM-00は開始条件を満たした。[開始審査](implementation-start-review.md)に対象・依存・未完了を記録。')
            s=s.replace('現時点では計算本体のコードを追加しない。','開始審査に合格する前は計算本体のコードを追加しない。現在IM-00は合格し、その範囲から着手可能。')
            s=s.replace('2026-11-03時点。草案・限定検証・採用完了を区別する。','以下の棚卸しは2026-11-03時点の履歴。最新の採用状態は冒頭と[開始審査](implementation-start-review.md)を優先する。草案・限定検証・採用完了を区別する。')
        elif name=='docs/requirements-evidence-map.md':
            a=s.find('E-SHELLの現在の追加証拠（2026-11-05）：')
            if a>=0:
                b=s.index('\n',a);s=s[:a]+'最新の計画は[A5-01〜07](practical-acceptance-plan.md)、[19ケース・37製品要求の台帳](verification/practical-acceptance-case-catalog.json)。1〜4の採用設計と最初の実装開始審査は完了。以下のE群は従来の追跡IDとして保持し、製品の未実装・未受入と方式採用を区別する。'+s[b:]
            lines=s.splitlines()
            groups={'E-AERO':'V-AERO','E-MULTI':'V-MULTI','E-SHELL':'V-PG','E-LAM':'V-MAT','E-PRE':'V-MAT／V-J1','E-JOINT':'V-JOINT','E-BEAM':'V-BEAM','E-SUPPORT':'V-LOAD','E-MAP':'V-MULTIの写像','E-STATIC':'X-S10-STATIC／P-CURVED-BOX','E-DYNAMIC':'V-J1／X-S10-DECAY','E-FLUTTER':'V-SEARCH／X-S10-FLUTTER','E-IMPORT':'V-USE','E-RESTART':'V-STATE','E-USE':'V-USE','E-PERF':'PERF-FL-01'}
            for i,line in enumerate(lines):
                for key,value in groups.items():
                    if line.startswith('| '+key+' |'):
                        cells=line.split('|');cells[-2]=' 最新の入力・参照・閾値・不足時の扱いはA5の'+value+'。採用設計は確定、製品移植・実務受入は未完了。S10実験認定と製品性能は未認定 ';lines[i]='|'.join(cells);break
            s='\n'.join(lines)+'\n'
        elif name=='docs/design-decisions.md':
            s=s.replace('具体的な理論・数値・I/Fの決定は未基準化であり、製品要求の合意と区別する。','最新の採用式・具体規約は構造／UVLM／J1／N4／IO4、受入計画はA5に基準化。最初のIM-00開始審査は合格。以下の初期棚卸し・追記はその経緯として保持する。')
            s=s.replace('具体数値解法／抽出器・データ実装は4、実務受入は5','具体数値解法／抽出器・型／I/FはN4／IO4で採用済み。受入計画はA5、製品移植・受入結果はIM-00〜08')
            lines=s.splitlines()
            for i,line in enumerate(lines):
                if re.match(r'^\| (NU|IF)-\d+ \|',line):
                    cells=line.split('|');cells[-2]=' N4-01〜10／IO4-01〜07で型・所有権・寿命・疎密解法・精度・失敗・再開・容量政策を採用。両OSの限定試作を監査。製品移植・全実規模性能・実務受入は未完了（[開始審査](implementation-start-review.md)） ';lines[i]='|'.join(cells)
            s='\n'.join(lines)+'\n'
        p.write_text(s,encoding='utf-8')
    for name,h in audit['sources'].items():assert sha(ROOT/name)==h,name
    decision=dict(passed=True,preimplementation_initial_package_ready=True,practical_acceptance_plan_complete=True,initial_implementation_review_complete=True,implementation_start_allowed_packages=['IM-00'],product_implementation_started=False,experimental_reference_qualified=False,product_performance_qualified=False,all_product_acceptance_complete=False,current_work='IM-00 product numerical representations, solvers and state contracts; not started',risks=[f'R-{i:02}' for i in range(1,8)])
    dp=ROOT/'docs/verification/implementation-readiness-disposition.json';assert not dp.exists();dp.write_text(json.dumps(decision,indent=2),encoding='utf-8')
    closure=dict(**decision,immutable={p.relative_to(ROOT).as_posix():sha(p) for p in [audit_path,dp,ROOT/'docs/implementation-start-review.md',ROOT/'docs/practical-acceptance-plan.md',ROOT/'docs/verification/practical-acceptance-plan-completion.md']},status_sync_source={Path(__file__).relative_to(ROOT).as_posix():sha(Path(__file__))},mutable_status_entries=fronts)
    cp=ROOT/'docs/verification/implementation-readiness-closure.json';assert not cp.exists();cp.write_text(json.dumps(closure,indent=2),encoding='utf-8');print(json.dumps(decision))
if __name__=='__main__':main()
