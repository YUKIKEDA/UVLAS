from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def main():
    fronts=['README.md','docs/preimplementation-plan.md','docs/structural-model-adoption-review.md','docs/requirements-evidence-map.md','docs/design-decisions.md','docs/roadmap.md']
    for name in fronts:
        p=ROOT/name;s=p.read_text(encoding='utf-8');base='docs/' if name=='README.md' else ''
        if name=='docs/roadmap.md':
            lines=s.splitlines()
            for i,line in enumerate(lines):
                if line.startswith('| 2.2 後流・組立・状態・荷重：採用完了 |'):lines[i]='| 2.2 後流・組立・状態・荷重：採用完了 | **16群が合格**。全相互組立・陰的先頭・自由／指定移流・全荷重とPG写像／接線・状態／再開を採用 | [仕様](uvlm-wake-adopted-specification.md)、[採用判定](verification/uvlm-wake-adoption-results.md)。製品の時間方式は3へ渡す |'
                if line.startswith('| **2.3 独立精度参照・誤差予算：現在地** |'):lines[i]='| **2.3 独立精度参照・誤差予算：現在地** | 定常2翼の独立総量比較と感度・空間指標は成立。非定常は独立時間／周波数と空間差を分離し追加中 | [実数値と残条件](verification/uvlm-independent-reference-progress.md)。非定常・帯の連続極限・誘導抗力・多面の独立比較を閉じる。AVL系列は再開しない |'
            s='\n'.join(lines)+'\n'
        s+=f'\n2026-11-05・2.3更新：[独立参照の実数値・未達・次の終了条件]({base}verification/uvlm-independent-reference-progress.md)。定常の総量・感度は限定成立、非定常の空間収束を追加中。UVLM全体と本体実装開始は未完了。\n'
        p.write_text(s,encoding='utf-8')
if __name__=='__main__':main()
