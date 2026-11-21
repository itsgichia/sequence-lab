"""Sequence parsing, statistics, translation and open reading frames."""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass

MAX_SEQ = 5000  # residues per sequence

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
DNA_LETTERS = set("ACGTN")
PROTEIN_LETTERS = set("ACDEFGHIKLMNPQRSTVWY*X")

# Average residue masses (Da), water added once per chain.
AA_MASS = {
    "A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.14, "E": 129.12, "Q": 128.13, "G": 57.05,
    "H": 137.14, "I": 113.16, "L": 113.16, "K": 128.17, "M": 131.19, "F": 147.18, "P": 97.12, "S": 87.08,
    "T": 101.10, "W": 186.21, "Y": 163.18, "V": 99.13,
}


class BioError(ValueError):
    pass


# ---------------------------------------------------------------- parsing
@dataclass
class Record:
    name: str
    seq: str


def clean(raw: str) -> str:
    """Strip whitespace, digits and a leading FASTA header; upper-case; U -> T."""
    lines = [ln for ln in raw.strip().splitlines() if not ln.startswith(">")]
    seq = re.sub(r"[\s\d]", "", "".join(lines)).upper().replace("U", "T")
    if not seq:
        raise BioError("No sequence found.")
    if len(seq) > MAX_SEQ:
        raise BioError(f"Sequence is longer than {MAX_SEQ:,} residues.")
    return seq


def parse_fasta(raw: str) -> list[Record]:
    records: list[Record] = []
    name, chunks = None, []
    for ln in raw.strip().splitlines():
        if ln.startswith(">"):
            if name is not None:
                records.append(Record(name, "".join(chunks)))
            name, chunks = ln[1:].strip() or f"seq{len(records) + 1}", []
        else:
            chunks.append(ln)
    if name is not None:
        records.append(Record(name, "".join(chunks)))
    if not records:
        raise BioError("No FASTA records found. Each sequence needs a '>' header line.")
    out = []
    for r in records:
        seq = re.sub(r"[\s\d]", "", r.seq).upper().replace("U", "T")
        if not seq:
            raise BioError(f"'{r.name}' has no sequence.")
        out.append(Record(r.name[:40], seq))
    return out


def kind_of(seq: str) -> str:
    letters = set(seq)
    if letters <= DNA_LETTERS:
        return "dna"
    if letters <= PROTEIN_LETTERS:
        return "protein"
    bad = sorted(letters - PROTEIN_LETTERS)
    raise BioError(f"Unexpected characters: {''.join(bad)[:10]}")


# ------------------------------------------------------------- sequence
def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq: str) -> str:
    return "".join(CODON_TABLE.get(seq[i : i + 3], "X") for i in range(0, len(seq) - 2, 3))


def find_orfs(seq: str, min_aa: int = 30) -> list[dict]:
    """ATG-to-stop open reading frames on both strands, longest first."""
    orfs = []
    n = len(seq)
    for strand, s in (("+", seq), ("-", reverse_complement(seq))):
        for frame in range(3):
            i = frame
            while i + 3 <= n:
                if s[i : i + 3] == "ATG":
                    j = i
                    while j + 3 <= n and CODON_TABLE.get(s[j : j + 3], "X") != "*":
                        j += 3
                    aa = (j - i) // 3
                    if aa >= min_aa:
                        start, end = (i, j + 3) if strand == "+" else (n - (j + 3), n - i)
                        orfs.append(
                            {
                                "strand": strand,
                                "frame": frame + 1,
                                "start": start + 1,  # 1-based, inclusive
                                "end": min(end, n),
                                "length_aa": aa,
                                "protein": translate(s[i : j + 3]).rstrip("*"),
                            }
                        )
                    i = j + 3
                else:
                    i += 3
    orfs.sort(key=lambda o: -o["length_aa"])
    return orfs[:20]


def stats(raw: str, min_orf_aa: int = 30) -> dict:
    seq = clean(raw)
    kind = kind_of(seq)
    counts = Counter(seq)
    n = len(seq)
    out = {"kind": kind, "length": n, "composition": [{"residue": k, "count": v, "percent": round(100 * v / n, 1)} for k, v in sorted(counts.items())]}
    if kind == "dna":
        gc = counts["G"] + counts["C"]
        at = counts["A"] + counts["T"]
        out["gc_percent"] = round(100 * gc / max(1, gc + at), 1)
        window = max(20, min(200, n // 20))
        out["gc_window"] = window
        out["gc_track"] = [
            {"pos": i + window // 2, "gc": round(100 * sum(1 for c in seq[i : i + window] if c in "GC") / window, 1)}
            for i in range(0, max(1, n - window + 1), max(1, window // 4))
        ]
        out["reverse_complement"] = reverse_complement(seq) if n <= 400 else None
        out["orfs"] = find_orfs(seq, min_orf_aa)
        out["frames"] = [{"frame": f + 1, "protein": translate(seq[f:])} for f in range(3)] if n <= 1200 else []
    else:
        mass = sum(AA_MASS.get(c, 110.0) for c in seq if c != "*") + 18.02
        out["mass_kda"] = round(mass / 1000, 2)
        hydrophobic = sum(counts[c] for c in "AVILMFWY")
        out["hydrophobic_percent"] = round(100 * hydrophobic / n, 1)
        out["charged_percent"] = round(100 * sum(counts[c] for c in "DEKRH") / n, 1)
    return out
