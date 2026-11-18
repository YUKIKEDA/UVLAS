"""Review regression: actual two-pass tree traversal and endpoint-based edge alias."""
import argparse
import json
from pathlib import Path


KINDS = ("cmp", "decode", "push", "pop", "visit")


def tree(boxes):
    nodes = []
    def build(items):
        bounds = tuple((min(boxes[i][k][0] for i in items), max(boxes[i][k][1] for i in items)) for k in range(3))
        index = len(nodes)
        nodes.append({"bounds": bounds, "count": len(items), "leaf": items[0] if len(items) == 1 else None})
        if len(items) > 1:
            mid = len(items) // 2
            nodes[index]["children"] = (build(items[:mid]), build(items[mid:]))
        return index
    build(list(range(len(boxes))))
    return nodes


def traversal(nodes):
    counts = dict.fromkeys(KINDS, 0)
    stack = [(0, 0)]
    counts["push"] = 1
    pairs = []
    while stack:
        i, j = stack.pop()
        counts["pop"] += 1
        counts["visit"] += 1
        a, b = nodes[i], nodes[j]
        if i == j and a["leaf"] is not None:
            continue
        overlap = True
        for k in range(3):
            for lo, hi in [(a["bounds"][k][0], b["bounds"][k][1]), (b["bounds"][k][0], a["bounds"][k][1])]:
                counts["cmp"] += 1
                counts["decode"] += 4
                if lo > hi:
                    overlap = False
                    break
            if not overlap:
                break
        if not overlap:
            continue
        if i == j:
            left, right = a["children"]
            children = [(left, left), (left, right), (right, right)]
        elif a["leaf"] is not None and b["leaf"] is not None:
            pairs.append(tuple(sorted((a["leaf"], b["leaf"]))))
            continue
        elif b["leaf"] is not None or (a["leaf"] is None and a["count"] >= b["count"]):
            children = [(c, j) for c in a["children"]]
        else:
            children = [(i, c) for c in b["children"]]
        stack.extend(reversed(children))
        counts["push"] += len(children)
    assert len(pairs) == len(set(pairs))
    return pairs, counts


def run_budget(nodes, limit, omit_second=False):
    n = len(nodes)
    v = n * (n + 1) // 2
    bounds = dict(zip(KINDS, (6 * v, 24 * v, 3 * v, v, v)))
    b = sum(bounds.values())  # All symbolic primitive weights are 1 in this witness.
    reserved = b if omit_second else 2 * b
    if reserved > limit:
        return {"success": False, "passes": 0, "reserved": reserved}
    first, count1 = traversal(nodes)
    assert all(count1[k] <= bounds[k] for k in KINDS)
    reserved += 7 * len(first) + 11  # Mock narrow/rest units, not product prices.
    if reserved > limit:
        return {"success": False, "passes": 1, "reserved": reserved}
    second, count2 = traversal(nodes)
    assert first == second and count1 == count2
    assert all(count2[k] <= bounds[k] for k in KINDS)
    return {"success": True, "passes": 2, "reserved": reserved, "candidates": len(first),
            "counts_per_pass": count1, "traversal_reserved": b if omit_second else 2 * b,
            "required_traversal_reservation": 2 * b}


def resolved_incidence(points, edges, edge_id, traversal_endpoints, omit_basis=False):
    owner_id = edges[edge_id]["owner"]
    canonical = tuple(points[p] for p in edges[owner_id]["endpoints"])
    raw = tuple(points[p] for p in edges[edge_id]["endpoints"])
    walked = tuple(points[p] for p in traversal_endpoints)
    assert canonical[0] < canonical[1]
    assert sorted(raw) == list(canonical) and sorted(walked) == list(canonical)
    raw_incidence = 1 if walked == raw else -1
    basis = 1 if raw == canonical else -1
    return raw_incidence if omit_basis else raw_incidence * basis


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Refusing to overwrite existing result")
    records = []
    for name, boxes, expected_k in [
        ("zero-candidates", [((0, 1), (0, 1), (0, 1)), ((3, 4), (3, 4), (3, 4))], 0),
        ("six-candidates", [((0, 1), (0, 1), (0, 1))] * 4, 6),
    ]:
        nodes = tree(boxes)
        v = len(nodes) * (len(nodes) + 1) // 2
        exact = 2 * 35 * v + expected_k * 7 + 11
        equal = run_budget(nodes, exact)
        short = run_budget(nodes, exact - 1)
        assert equal["success"] and equal["passes"] == 2 and equal["candidates"] == expected_k
        assert not short["success"] and short["passes"] == 1
        mutant = run_budget(nodes, exact, omit_second=True)
        assert mutant["success"]
        assert mutant["traversal_reserved"] != mutant["required_traversal_reservation"]
        # Missing second reservation also wrongly admits the one-short input.
        assert run_budget(nodes, exact - 1, omit_second=True)["success"]
        records.append({"case": name, "equal": equal, "one_short": short, "missing_second_detected": True})
    points = {2: 2, 3: 3, 20: 2, 30: 3}
    edges = {100: {"owner": 100, "endpoints": (2, 3)}, 200: {"owner": 100, "endpoints": (30, 20)}}
    a = resolved_incidence(points, edges, 200, (20, 30))
    b = resolved_incidence(points, edges, 100, (3, 2))
    wake = -resolved_incidence(points, edges, 200, (30, 20))
    t = -wake * (-a)
    assert [a, b, a + b, wake, t] == [1, -1, 0, 1, 1]
    wrong_a = resolved_incidence(points, edges, 200, (20, 30), omit_basis=True)
    wrong_t = -wake * (-wrong_a)
    assert [wrong_a, wrong_a + b, wrong_t] == [-1, -2, -1]
    assert [wrong_a, wrong_a + b, wrong_t] != [1, 0, 1]
    result = {"work": records, "point_owners": points, "edges": edges,
              "alias": {"geometric_incidence": [a, b], "sum": a + b, "wake_induction": wake, "T": t,
                        "without_basis": {"incidence": wrong_a, "sum": wrong_a + b, "T": wrong_t}},
              "tolerance": "0", "negative_controls": 3}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    print("two work cases / equal and one-short / reversed alias / 3 negative controls: passed")


if __name__ == "__main__":
    main()
