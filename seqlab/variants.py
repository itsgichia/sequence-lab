"""Call variants between a reference coding sequence and a sample, and predict their effect.

The reference is a coding sequence (CDS) read from the first base of the start
codon, so a position `p` is HGVS `c.p`. Variants are called from a global
alignment, then each one is applied to the reference on its own and both
coding sequences are translated to classify the consequence:

    synonymous, missense, nonsense (stop gained), stop lost, start lost,
    in-frame insertion or deletion, frameshift.

Notation follows HGVS: `c.76A>G`, `c.120_122del`, `c.45_46insT`, and
`p.Ile26Val`, `p.Trp40Ter`, `p.Lys15fs`, `p.Gly41del`.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

from .align import align
from .sequence import BioError, CODON_TABLE, clean, kind_of, translate

THREE_LETTER = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "E": "Glu", "Q": "Gln", "G": "Gly",
    "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
    "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter", "X": "Xaa",
}
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


@dataclass
class Variant:
    """One difference between reference and sample, in reference coordinates."""

    kind: str  # substitution | deletion | insertion | delins
    start: int  # 1-based reference position (for insertions: the base before the insertion)
    end: int  # 1-based inclusive
    ref: str
    alt: str
    cdna: str = ""
    protein: str = ""
    consequence: str = ""
    codon: int | None = None
    transition: bool | None = None

    def apply(self, reference: str) -> str:
        """Return the reference with only this variant applied."""
        if self.kind == "insertion":
            return reference[: self.start] + self.alt + reference[self.start :]
        return reference[: self.start - 1] + self.alt + reference[self.end :]


def _hgvs_c(v: Variant) -> str:
    if v.kind == "substitution":
        return f"c.{v.start}{v.ref}>{v.alt}"
    if v.kind == "insertion":
        return f"c.{v.start}_{v.start + 1}ins{v.alt}"
    span = f"{v.start}" if v.start == v.end else f"{v.start}_{v.end}"
    if v.kind == "deletion":
        return f"c.{span}del"
    return f"c.{span}delins{v.alt}"


def _aa(code: str) -> str:
    return THREE_LETTER.get(code, "Xaa")


def _consequence(v: Variant, reference: str) -> tuple[str, str]:
    """Classify a variant and give its protein-level HGVS description."""
    ref_prot = translate(reference)
    alt_prot = translate(v.apply(reference))
    ref_len = len(v.ref) if v.kind != "insertion" else 0
    alt_len = len(v.alt)
    frame_change = (alt_len - ref_len) % 3 != 0

    # First differing codon.
    first = next((k for k, (x, y) in enumerate(zip(ref_prot, alt_prot)) if x != y), None)
    if first is None:
        if len(alt_prot) == len(ref_prot):
            codon_index = (v.start - 1) // 3
            res = ref_prot[codon_index] if codon_index < len(ref_prot) else "X"
            return "synonymous", f"p.{_aa(res)}{codon_index + 1}="
        first = min(len(ref_prot), len(alt_prot))

    pos = first + 1
    ref_res = ref_prot[first] if first < len(ref_prot) else ""
    alt_res = alt_prot[first] if first < len(alt_prot) else ""

    if pos == 1 and ref_res == "M" and alt_res != "M":
        return "start lost", "p.Met1?"
    if frame_change:
        # Distance to the new stop codon, if one appears within the frame-shifted protein.
        stop = alt_prot.find("*", first)
        ter = f"Ter{stop - first + 1}" if stop >= 0 else "Ter?"
        return "frameshift", f"p.{_aa(ref_res or 'X')}{pos}{_aa(alt_res) if alt_res else ''}fs{ter}"
    if v.kind == "substitution":
        if alt_res == "*":
            return "nonsense", f"p.{_aa(ref_res)}{pos}Ter"
        if ref_res == "*":
            ext = alt_prot.find("*", first)
            return "stop lost", f"p.Ter{pos}{_aa(alt_res)}ext*{ext - first if ext >= 0 else '?'}"
        return "missense", f"p.{_aa(ref_res)}{pos}{_aa(alt_res)}"
    # In-frame indel: describe the run of residues that changed.
    n_ref = ref_len // 3 if v.kind != "insertion" else 0
    n_alt = alt_len // 3
    if v.kind == "deletion" or (v.kind == "delins" and n_alt < n_ref):
        lost = ref_prot[first : first + max(1, n_ref - n_alt)]
        tail = "" if len(lost) == 1 else f"_{_aa(lost[-1])}{first + len(lost)}"
        label = f"p.{_aa(lost[0])}{pos}{tail}del"
        if "*" in alt_prot[first : first + n_alt + 1] and "*" not in ref_prot[first : first + n_ref + 1]:
            return "in-frame deletion, stop gained", label
        return "in-frame deletion", label
    if v.kind == "insertion" or n_alt > n_ref:
        gained = alt_prot[first : first + n_alt - n_ref]
        if "*" in gained:
            return "in-frame insertion, stop gained", f"p.{_aa(ref_res)}{pos}Ter"
        left = ref_prot[first - 1] if first > 0 else ""
        inserted = "".join(_aa(c) for c in gained)
        return "in-frame insertion", f"p.{_aa(left)}{first}_{_aa(ref_res)}{pos}ins{inserted}"
    return "in-frame substitution", f"p.{_aa(ref_res)}{pos}{_aa(alt_res)}"


def call_variants(reference_raw: str, sample_raw: str, match: int = 2, mismatch: int = -1, gap: int = -3) -> dict:
    """Align a sample coding sequence to its reference and describe every difference."""
    reference, sample = clean(reference_raw), clean(sample_raw)
    if kind_of(reference) != "dna" or kind_of(sample) != "dna":
        raise BioError("Variant calling needs two DNA sequences.")
    if not reference.startswith("ATG"):
        raise BioError("The reference should be a coding sequence starting with ATG (c.1 is the A of the start codon).")
    al = align(reference, sample, "global", match=match, mismatch=mismatch, gap=gap)

    variants: list[Variant] = []
    ref_pos = 0  # number of reference bases consumed
    run: Variant | None = None

    def flush():
        nonlocal run
        if run is not None:
            variants.append(run)
            run = None

    for x, y in zip(al["a"], al["b"]):
        if x != "-":
            ref_pos += 1
        if x == y:
            flush()
            continue
        if x == "-":  # insertion in the sample after ref_pos
            if run and run.kind in ("insertion", "delins") and run.end == ref_pos:
                run.alt += y
            elif run and run.kind == "deletion" and run.end == ref_pos:
                run.kind, run.alt = "delins", run.alt + y
            else:
                flush()
                run = Variant("insertion", ref_pos, ref_pos, "", y)
            continue
        if y == "-":  # deletion of reference base ref_pos
            if run and run.kind == "deletion" and run.end == ref_pos - 1:
                run.end, run.ref = ref_pos, run.ref + x
            elif run and run.kind in ("substitution", "delins") and run.end == ref_pos - 1:
                run.kind, run.end, run.ref = "delins", ref_pos, run.ref + x
            else:
                flush()
                run = Variant("deletion", ref_pos, ref_pos, x, "")
            continue
        # substitution
        flush()
        run = Variant("substitution", ref_pos, ref_pos, x, y)
    flush()

    for v in variants:
        v.cdna = _hgvs_c(v)
        v.consequence, v.protein = _consequence(v, reference)
        v.codon = (v.start - 1) // 3 + 1 if v.kind != "insertion" else v.start // 3 + 1
        if v.kind == "substitution":
            v.transition = (v.ref, v.alt) in TRANSITIONS

    subs = [v for v in variants if v.kind == "substitution"]
    ti = sum(1 for v in subs if v.transition)
    tv = len(subs) - ti
    ref_prot, sam_prot = translate(reference), translate(sample)
    counts: dict[str, int] = {}
    for v in variants:
        counts[v.consequence] = counts.get(v.consequence, 0) + 1
    return {
        "reference_length": len(reference),
        "sample_length": len(sample),
        "identity_percent": al["identity_percent"],
        "variants": [asdict(v) for v in variants],
        "n_variants": len(variants),
        "substitutions": len(subs),
        "transitions": ti,
        "transversions": tv,
        "ti_tv": round(ti / tv, 2) if tv else None,
        "indels": len(variants) - len(subs),
        "consequences": counts,
        "reference_protein": ref_prot,
        "sample_protein": sam_prot,
        "protein_truncated": sam_prot.find("*") != -1 and sam_prot.find("*") < max(0, ref_prot.find("*")),
    }


__all__ = ["Variant", "call_variants", "CODON_TABLE"]
