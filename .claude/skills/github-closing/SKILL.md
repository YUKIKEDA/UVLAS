---
name: github-closing
description: "PRが終えていないIssueを、開いたままにする。PRの本文、commitのメッセージ、gh pr create、Closes、Issue番号を名指しする文を書くときに使う。"
---

# 閉じる語

GitHubは、PRの題、本文、squashのcommitメッセージの中で、閉じる動詞がIssue番号のすぐ前にあると、merge時にそのIssueを完了にする。否定でも一致する：`This does not close #12` は #12 を完了にする。

動詞：`close`、`closes`、`closed`、`fix`、`fixes`、`fixed`、`resolve`、`resolves`、`resolved`。

## Issueを終える

`## Related` の下に、一行、他に何も書かずに：

```text
Closes #12
```

同じ番号の隣で、その動詞を他の場所に繰り返さない。commitの本文には決して書かない。

## Issueを開いたままにする

`#12` を、それらの動詞を前に置かずに書く。

```text
Part of #12. 完了gateの「両OSのクリーンビルド」はまだ満たしていない。
```

親Issue（実装パッケージ）を指す子のPRは、常に `Part of` である。

## 送る前に

1. 題と本文の全体を読む。Issue番号ごとに、このPRがそれを終えるかを決める。Issueの完了条件の各行を、行ったことの横に置いて決める。
2. 終えないのに、番号の前の語が閉じる動詞なら、文を書き直す。
3. `gh pr create`（または `gh pr edit`）を走らせる。
4. `gh pr view <n> --json closingIssuesReferences` を走らせ、一覧を手順1と比べる。

別のブランチに積んだPRは、閉じる参照を示さない。GitHubは、既定のブランチへ向くPRにだけ結び付ける。PRが `main` を向いたときに、一覧をもう一度確かめる。
