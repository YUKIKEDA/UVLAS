"""Synchronize the canonical roadmap after, and only after, candidate qualification."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[2]

def main():
    assert json.loads((ROOT/'docs/verification/paired-geometry-ancf-candidate-result.json').read_text())['adoption_candidate_qualified']
    p=ROOT/'docs/roadmap.md';lines=p.read_text(encoding='utf-8-sig').splitlines()
    replacements={
    '**現在地は':'**現在地は1.1主シェル。PG-ANCFを最終採用審査へ進める候補として選定し、R2の具体案とR3の今回の適格性検査を終了した。次はR4の一般面・物理PDE・実務力学との整合、その後R5の主シェル最終採用判断。** [候補判定・主・補足16配置48右辺](verification/paired-geometry-ancf-candidate-results.md)、[数式と処理境界](paired-geometry-ancf-definition.md)。旧MMS-1と独立r/d補完GCの不採用は維持。主シェル最終採用・本体実装開始は未承認。',
    '現MMS-1の不採用は':'現MMS-1の不採用は取り消さない。PG-ANCFでは位置と厚さ勾配を一緒に持ち上げるスカラー形状関数に変え、同じ全3Dエネルギーと基準で選別した。ANCF要件、実務範囲、1%・薄肉増幅2倍を含む既定基準を維持する。',
    '| **R2 再定式化':'| **R2 具体案：終了** | [PG-ANCF](paired-geometry-ancf-definition.md)：φ₁=−κC−zI₀、φ₃=b−κxzI₀を各3Cartesian方向に持つ対の位置補完。基本r/d/cは独立、線形S・一定質量・全3Dエネルギーを維持。一般面への参照Ritz手順と主要処理境界も定義 | 独立r/d補完の未達を変える具体式を得た。一般面の手順を定義したことと、一般面での品質成立は区別する。後者はR4 |',
    '| R3 決定的':'| **R3 決定的有限厚さ検査：今回の範囲を終了** | [正則／歪み・2厚さ・標準8配置＋独立曲率8配置](verification/paired-geometry-ancf-candidate-results.md)、主・補足16配置48右辺、別観測・別SVD・原停留・共有・質量・剛体／仕事／接線を通過 | 1%・薄肉増幅2倍を同時に満たし、採用審査候補として選定。数値床の標準場だけで精度を主張しない。物理PDE・一般面・材料等の最終採用は未完了 |',
    '| R4 実務':'| **R4 実務に必要な力学との整合：次** | 一般面の参照Ritz持上げと必要二階微分／全トレース、独立な物理体積／面荷重・反力・仕事、採用材料／積層・初期ひずみ、非線形U/f/K、慣性／固有特性・再開、代表規模の容量をPGの同じSへ接続する | 採用を左右する代表照合と適用／拒否条件、具体I/F、実装後へ残す検証の入力・参照・許容差を対応させる。一般面と費用の未達なら候補を最終採用しない |',
    'R1〜R5は':'R1〜R5は実装前6作業の「1」の中で主シェルを閉じる作業。R2/R3で初めて今回の判定範囲を満たすPG候補に絞れた。残るR4/R5の仕事量は同じでないため「3/5=60%」とは換算しない。次は上表のR4を一組で進め、一般面・物理PDE・材料／慣性と容量の採用リスクを閉じる。採用日を保証できる段階ではない。',
    '### R1の実施範囲':'### R1の実施範囲とR2への引継ぎ（終了した設計入力）',
    'R2の設計対象は':'R1が要求したのは、同じ離散曲げ場が厚さ変動の曲げ近似、平均応力仕事、鉛直面共有、上下面整合を同時に満たす持上げ。これに対してPGの対の位置補完を定義し、R3で1%・薄肉増幅2倍と数値条件を別に判定した。[今回の採否](verification/paired-geometry-ancf-candidate-results.md)。以下はR1実施前の固定範囲を保持する履歴であり、現在の作業指示ではない。',
    '- 次の試作はR2':'- R2/R3で今回の候補選別を閉じた。次はR4の一般面・物理PDE・材料／慣性・容量を一組で判定し、候補の局所監査や既に通過した条件を最初からやり直さない。',
    '| **1 構造モデル採用**':'| **1 構造モデル採用** | **進行中・主シェル未採用**。PG-ANCFを採用審査候補に選定しR2/R3の今回の範囲を終了 | R4の一般面・物理PDE・実務力学／容量→R5採用。その後、梁・慣性・接続等を統合 | ST-01〜10/CV-03の採用式・適用条件・必要微分・接続・独立参照と検証計画が揃う |',
    '| 1.1 主シェル':'| 1.1 主シェルの補間・ロッキング・求積・幾何 | **PG-ANCFを採用審査候補に選定**。R2具体案、R3の主・補足16配置・48右辺と独立照合を完了 | R4一般面・物理PDE・実務力学／容量と、R5の最終採用は未完了。[候補判定](verification/paired-geometry-ancf-candidate-results.md) |',
    '以下は各実施時点の履歴。':'以下は各実施時点の履歴。履歴中の「次」は当時の予定であり、現在の作業指示は冒頭の「いま、どこにいるか」と「主シェルを閉じるまで」を正とする。',
    '| 次の採用判断':'| 次の採用判断 | PG候補の一般曲面・物理PDE・材料／積層・非線形・慣性・費用をR4で統合し、R5で最終採用案を出す | 候補の成立と主シェル最終採用を区別し、未達があれば基準を緩めず不採用とする |',
    }
    found={k:0 for k in replacements}
    for i,line in enumerate(lines):
        for prefix,text in replacements.items():
            if line.startswith(prefix):lines[i]=text;found[prefix]+=1;break
    assert all(v==1 for v in found.values()),found
    p.write_text('\n'.join(lines)+'\n',encoding='utf-8');print('Roadmap synchronized:',len(found),'entries')

if __name__=='__main__':main()
