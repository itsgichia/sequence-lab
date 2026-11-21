"""Synthetic example data with known answers, so every tool can be checked by eye."""

from __future__ import annotations

import random

from .align import MAX_CELLS
from .sequence import CODON_TABLE, MAX_SEQ, parse_fasta, translate
from .tree import MAX_TREE_LEN, MAX_TREE_SEQS


def example_family(seed: int = 7) -> str:
    """A synthetic gene family: one ancestor, two clades, point mutations and one deletion.

    Synthetic so the "true" tree is known: (A1, A2, A3) and (B1, B2, B3) are the two clades,
    with C an outgroup. Useful for seeing what neighbour joining recovers.
    """
    rng = random.Random(seed)
    ancestor = "ATG" + "".join(rng.choice("ACGT") for _ in range(297))

    def mutate(seq: str, k: int) -> str:
        s = list(seq)
        for pos in rng.sample(range(3, len(s)), k):
            s[pos] = rng.choice([c for c in "ACGT" if c != s[pos]])
        return "".join(s)

    clade_a, clade_b = mutate(ancestor, 18), mutate(ancestor, 18)
    seqs = {
        "A1 demo": mutate(clade_a, 4),
        "A2 demo": mutate(clade_a, 5),
        "A3 demo": mutate(clade_a, 9),
        "B1 demo": mutate(clade_b, 4),
        "B2 demo": mutate(clade_b, 6),
        "B3 demo": mutate(clade_b, 6)[:150] + mutate(clade_b, 6)[162:],  # a 12 bp deletion
        "C outgroup": mutate(ancestor, 45),
    }
    return "\n".join(f">{k}\n{v}" for k, v in seqs.items())


def example_gene(seed: int = 11) -> str:
    """A synthetic 240 bp coding sequence with an open frame and a final stop codon."""
    rng = random.Random(seed)
    codons = [c for c, aa in CODON_TABLE.items() if aa != "*"]
    body = "".join(rng.choice(codons) for _ in range(78))
    return "ATG" + body + "TAA"


def example_variants() -> tuple[str, str]:
    """Reference and sample coding sequences differing by one variant of each class.

    Engineered on the synthetic gene so the expected calls are known:
    a synonymous change, a missense change, a nonsense change, an in-frame
    3 bp deletion and a 1 bp insertion that shifts the frame near the end.
    """
    ref = example_gene()
    protein = translate(ref)
    s = list(ref)

    def codon_at(i: int) -> str:
        return "".join(s[3 * i : 3 * i + 3])

    # Synonymous: codon 5, change the third base to another codon for the same residue.
    i = 4
    for alt in ("A", "C", "G", "T"):
        trial = codon_at(i)[:2] + alt
        if trial != codon_at(i) and CODON_TABLE[trial] == protein[i]:
            s[3 * i + 2] = alt
            break
    # Missense: codon 12, first base to a different non-stop residue.
    i = 11
    for alt in ("A", "C", "G", "T"):
        trial = alt + codon_at(i)[1:]
        if trial != codon_at(i) and CODON_TABLE[trial] not in ("*", protein[i]):
            s[3 * i] = alt
            break
    # Nonsense: codon 60 -> TAG (do this late so the earlier calls are still in-frame).
    i = 59
    s[3 * i : 3 * i + 3] = list("TAG")
    # In-frame deletion of codon 25 (bases 73-75).
    i = 24
    del s[3 * i : 3 * i + 3]
    # Frameshift: insert a T after base 200 of the (now 3 bp shorter) sequence, reference coordinate 203.
    s.insert(203, "T")
    return ref, "".join(s)


def examples() -> dict:
    fam = example_family()
    records = parse_fasta(fam)
    ref, sample = example_variants()
    return {
        "sequence": f">A1 demo (synthetic 300 bp gene)\n{records[0].seq}",
        "align_a": f">A1 demo\n{records[0].seq}",
        "align_b": f">B3 demo\n{records[5].seq}",
        "fasta": fam,
        "variant_reference": f">reference CDS (synthetic, 240 bp)\n{ref}",
        "variant_sample": f">sample CDS\n{sample}",
        "pgx": {"CYP2C19": ["*1", "*2"], "CYP2D6": ["*4", "*10"], "TPMT": ["*1", "*3A"], "SLCO1B1": ["T", "C"]},
        "limits": {
            "max_seq": MAX_SEQ,
            "max_tree_seqs": MAX_TREE_SEQS,
            "max_tree_len": MAX_TREE_LEN,
            "max_align_cells": MAX_CELLS,
        },
    }
