"""Command line interface.

    seqlab stats gene.fa
    seqlab align a.fa b.fa --mode local
    seqlab tree family.fa --distance jc
    seqlab variants reference.fa sample.fa
    seqlab pgx CYP2C19 '*1' '*2'
    seqlab examples --out examples/      # writes the bundled example files
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import align, call_variants, examples, report, stats, tree
from .sequence import BioError


def _read(path: str) -> str:
    return sys.stdin.read() if path == "-" else Path(path).read_text()


def _emit(obj: dict, as_json: bool) -> None:
    if as_json:
        json.dump(obj, sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        for k, v in obj.items():
            if isinstance(v, (list, dict)):
                v = json.dumps(v)[:200] + ("..." if len(json.dumps(v)) > 200 else "")
            print(f"{k:>20}: {v}")


def _print_alignment(al: dict, width: int = 60) -> None:
    for i in range(0, al["length"], width):
        print(f"a  {al['a'][i:i + width]}\n   {al['midline'][i:i + width]}\nb  {al['b'][i:i + width]}\n")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="seqlab", description=__doc__.split("\n\n")[0])
    p.add_argument("--json", action="store_true", help="print full JSON output")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("stats", help="composition, GC content, ORFs and translation")
    s.add_argument("file")
    s.add_argument("--min-orf", type=int, default=30)

    a = sub.add_parser("align", help="pairwise alignment")
    a.add_argument("a")
    a.add_argument("b")
    a.add_argument("--mode", choices=("global", "local"), default="global")
    a.add_argument("--match", type=int, default=2)
    a.add_argument("--mismatch", type=int, default=-1)
    a.add_argument("--gap", type=int, default=-2)

    t = sub.add_parser("tree", help="neighbour-joining tree from a FASTA file")
    t.add_argument("file")
    t.add_argument("--distance", choices=("jc", "p"), default="jc")

    v = sub.add_parser("variants", help="call variants between a reference CDS and a sample")
    v.add_argument("reference")
    v.add_argument("sample")

    g = sub.add_parser("pgx", help="phenotype and drug recommendations for a diplotype")
    g.add_argument("gene")
    g.add_argument("allele_1")
    g.add_argument("allele_2")

    e = sub.add_parser("examples", help="write the bundled example files")
    e.add_argument("--out", default="examples")

    args = p.parse_args(argv)
    try:
        if args.cmd == "stats":
            _emit(stats(_read(args.file), args.min_orf), args.json)
        elif args.cmd == "align":
            al = align(_read(args.a), _read(args.b), args.mode, args.match, args.mismatch, args.gap)
            if args.json:
                _emit(al, True)
            else:
                print(f"{al['mode']} alignment, score {al['score']}, identity {al['identity_percent']}%, gaps {al['gaps']}\n")
                _print_alignment(al)
        elif args.cmd == "tree":
            out = tree(_read(args.file), args.distance)
            _emit(out, True) if args.json else print(out["newick"])
        elif args.cmd == "variants":
            out = call_variants(_read(args.reference), _read(args.sample))
            if args.json:
                _emit(out, True)
            else:
                print(f"{out['n_variants']} variants, Ti/Tv {out['ti_tv']}, consequences {out['consequences']}\n")
                for var in out["variants"]:
                    print(f"{var['cdna']:<18} {var['protein']:<24} {var['consequence']}")
        elif args.cmd == "pgx":
            out = report(args.gene, args.allele_1, args.allele_2)
            if args.json:
                _emit(out, True)
            else:
                print(f"{out['gene']} {out['diplotype']}: {out['phenotype']}")
                if out["activity_score"] is not None:
                    print(f"activity score {out['activity_score']}")
                for r in out["recommendations"]:
                    print(f"\n{r['drug']}: {r['recommendation']}\n  ({r['guideline']})")
        elif args.cmd == "examples":
            out = Path(args.out)
            out.mkdir(parents=True, exist_ok=True)
            ex = examples()
            (out / "gene.fa").write_text(ex["sequence"] + "\n")
            (out / "family.fa").write_text(ex["fasta"] + "\n")
            (out / "reference_cds.fa").write_text(ex["variant_reference"] + "\n")
            (out / "sample_cds.fa").write_text(ex["variant_sample"] + "\n")
            print(f"wrote 4 files to {out}/")
    except BioError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
