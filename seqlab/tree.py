"""Distance matrices, neighbour joining, Newick export and a plotting layout."""

from __future__ import annotations

import numpy as np

from .align import align
from .sequence import BioError, kind_of, parse_fasta

MAX_TREE_SEQS = 12
MAX_TREE_LEN = 1500


def p_distance(a: str, b: str) -> float:
    al = align(a, b, "global")
    pairs = [(x, y) for x, y in zip(al["a"], al["b"]) if x != "-" and y != "-"]
    if not pairs:
        return 1.0
    return sum(1 for x, y in pairs if x != y) / len(pairs)


def jukes_cantor(p: float) -> float:
    p = min(p, 0.7499)
    return -0.75 * float(np.log(1 - 4 * p / 3))


def neighbour_joining(names: list[str], D: np.ndarray) -> dict:
    """Saitou & Nei neighbour joining. Returns a nested {name, length, children} tree."""
    nodes = [{"name": nm, "children": []} for nm in names]
    D = D.astype(float).copy()
    active = list(range(len(names)))
    while len(active) > 2:
        k = len(active)
        sub = D[np.ix_(active, active)]
        r = sub.sum(axis=1)
        Q = (k - 2) * sub - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = np.unravel_index(int(Q.argmin()), Q.shape)
        d_ij = sub[i, j]
        li = 0.5 * d_ij + (r[i] - r[j]) / (2 * (k - 2))
        lj = d_ij - li
        a, b = active[i], active[j]
        nodes[a]["length"], nodes[b]["length"] = max(li, 0.0), max(lj, 0.0)
        new = {"name": "", "children": [nodes[a], nodes[b]]}
        nodes.append(new)
        u = len(nodes) - 1
        D = np.pad(D, ((0, 1), (0, 1)))
        for x in active:
            if x in (a, b):
                continue
            D[u, x] = D[x, u] = 0.5 * (D[a, x] + D[b, x] - d_ij)
        active = [x for x in active if x not in (a, b)] + [u]
    a, b = active
    d = max(D[a, b], 0.0)
    nodes[a]["length"], nodes[b]["length"] = d / 2, d / 2
    return {"name": "", "children": [nodes[a], nodes[b]], "length": 0.0}


def newick(node: dict) -> str:
    def rec(n):
        if n["children"]:
            inner = ",".join(rec(c) for c in n["children"])
            return f"({inner}):{n.get('length', 0):.4f}" if n.get("length") is not None else f"({inner})"
        return f"{n['name'].replace(' ', '_')}:{n.get('length', 0):.4f}"

    return rec(node).rsplit(":", 1)[0] + ";"


def layout(node: dict) -> list[dict]:
    """Flatten the tree into rows with x (distance from root) and y (leaf order) for drawing."""
    rows, y = [], [0]

    def rec(n, x, parent):
        idx = len(rows)
        rows.append({"id": idx, "parent": parent, "name": n["name"], "x": x, "y": None, "leaf": not n["children"]})
        if n["children"]:
            ys = [rec(c, x + c.get("length", 0.0), idx) for c in n["children"]]
            rows[idx]["y"] = sum(ys) / len(ys)
        else:
            rows[idx]["y"] = y[0]
            y[0] += 1
        return rows[idx]["y"]

    rec(node, 0.0, None)
    return rows


def tree(fasta: str, distance: str = "jc") -> dict:
    records = parse_fasta(fasta)
    if len(records) < 3:
        raise BioError("A tree needs at least three sequences.")
    if len(records) > MAX_TREE_SEQS:
        raise BioError(f"At most {MAX_TREE_SEQS} sequences, please.")
    for r in records:
        if len(r.seq) > MAX_TREE_LEN:
            raise BioError(f"'{r.name}' is longer than {MAX_TREE_LEN:,} residues.")
    kinds = {kind_of(r.seq) for r in records}
    if len(kinds) > 1:
        raise BioError("Mix of DNA and protein sequences.")
    kind = kinds.pop()
    n = len(records)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            p = p_distance(records[i].seq, records[j].seq)
            D[i, j] = D[j, i] = jukes_cantor(p) if (distance == "jc" and kind == "dna") else p
    root = neighbour_joining([r.name for r in records], D)
    return {
        "kind": kind,
        "n": n,
        "distance": "Jukes-Cantor" if (distance == "jc" and kind == "dna") else "p-distance",
        "names": [r.name for r in records],
        "matrix": [[round(float(v), 4) for v in row] for row in D],
        "newick": newick(root),
        "nodes": layout(root),
    }
