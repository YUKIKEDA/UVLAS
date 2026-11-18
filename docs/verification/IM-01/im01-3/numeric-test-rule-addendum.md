不変量の残差には、製品から得た値と微分を使う。2026-11-17のIssue #108では、当初の膜残差をoracleの一次微分だけで評価していた。製品の係数から得た一次微分を独立oracleへ直接照合し、その実際の微分で膜残差を検査するよう補強した。oracle自身の残差だけを製品の不変量検証と数えない。

```cpp
// BAD: only checking the oracle's own invariant
CHECK(membrane(reference_jet) == 0);
// GOOD: compare independently, then check the actual product output
CHECK(actual_jet == reference_jet);
CHECK(membrane(actual_jet) == 0);
```

