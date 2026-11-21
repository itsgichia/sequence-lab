import pytest

from seqlab import (
    BioError,
    align,
    call_variants,
    catalogue,
    example_variants,
    examples,
    panel,
    phenotype,
    report,
    reverse_complement,
    stats,
    translate,
    tree,
)
from seqlab.cli import main


# ---------------------------------------------------------------- sequence
def test_translate_and_reverse_complement():
    assert translate("ATGGCCTAA") == "MA*"
    assert reverse_complement("ATGC") == "GCAT"


def test_stats_on_example_dna():
    out = stats(examples()["sequence"])
    assert out["kind"] == "dna" and out["length"] == 300
    assert 0 < out["gc_percent"] < 100
    assert out["orfs"] and out["orfs"][0]["strand"] in "+-"
    assert len(out["frames"]) == 3


def test_stats_on_protein():
    out = stats("MKTAYIAKQRQISFVKSHFSRQ")
    assert out["kind"] == "protein" and out["mass_kda"] > 2


def test_bad_characters_rejected():
    with pytest.raises(BioError):
        stats("ACGT123!!")


# --------------------------------------------------------------- alignment
def test_needleman_wunsch_textbook_example():
    out = align("GATTACA", "GCATGCU", "global", match=1, mismatch=-1, gap=-1)
    assert out["score"] == 0
    assert out["length"] == 8


def test_smith_waterman_finds_shared_core():
    out = align("ACGTACGTTT", "TTACGTAC", "local", match=2, mismatch=-1, gap=-2)
    assert out["a"] == "ACGTAC" and out["b"] == "ACGTAC" and out["score"] == 12


def test_align_rejects_mixed_types():
    with pytest.raises(BioError):
        align("ACGT", "MKT")


# -------------------------------------------------------------------- tree
def test_tree_recovers_synthetic_clades():
    out = tree(examples()["fasta"])
    assert out["n"] == 7 and len(out["nodes"]) == 13
    names, D = out["names"], out["matrix"]
    idx = {n: i for i, n in enumerate(names)}
    within = D[idx["A1 demo"]][idx["A2 demo"]]
    between = D[idx["A1 demo"]][idx["B1 demo"]]
    assert within < between
    assert out["newick"].endswith(";") and "C_outgroup" in out["newick"]


def test_tree_needs_three():
    with pytest.raises(BioError):
        tree(">a\nACGT\n>b\nACGA")


# ---------------------------------------------------------------- variants
def test_variant_classes_are_called_correctly():
    ref, sample = example_variants()
    out = call_variants(ref, sample)
    cons = [v["consequence"] for v in out["variants"]]
    assert cons == ["synonymous", "missense", "in-frame deletion", "nonsense", "frameshift"]
    v = {c["consequence"]: c for c in out["variants"]}
    assert v["synonymous"]["protein"].endswith("=")
    assert v["nonsense"]["protein"].endswith("Ter")
    assert "fs" in v["frameshift"]["protein"]
    assert v["in-frame deletion"]["cdna"] == "c.73_75del"
    assert v["nonsense"]["cdna"].startswith("c.") and v["nonsense"]["codon"] == 60
    assert out["protein_truncated"]
    assert out["substitutions"] == 3 and out["indels"] == 2


def test_single_substitution_notation():
    ref = "ATGAAAGGGTAA"  # M K G *
    out = call_variants(ref, "ATGAAAGGCTAA")  # GGG -> GGC still Gly
    assert out["variants"][0]["cdna"] == "c.9G>C"
    assert out["variants"][0]["protein"] == "p.Gly3="
    out = call_variants(ref, "ATGTAAGGGTAA")  # AAA -> TAA stop
    assert out["variants"][0]["consequence"] == "nonsense" and out["variants"][0]["protein"] == "p.Lys2Ter"


def test_variants_need_cds():
    with pytest.raises(BioError):
        call_variants("CCCAAA", "CCCAAA")


# --------------------------------------------------------------------- pgx
def test_cyp2c19_phenotypes():
    assert phenotype("CYP2C19", "*1", "*1")["phenotype"] == "Normal metaboliser"
    assert phenotype("CYP2C19", "*2", "*2")["phenotype"] == "Poor metaboliser"
    assert phenotype("CYP2C19", "*1", "*2")["phenotype"] == "Intermediate metaboliser"
    assert phenotype("CYP2C19", "*2", "*17")["phenotype"] == "Intermediate metaboliser"
    assert phenotype("CYP2C19", "*1", "*17")["phenotype"] == "Rapid metaboliser"
    assert phenotype("CYP2C19", "*17", "*17")["phenotype"] == "Ultrarapid metaboliser"


def test_cyp2d6_activity_score():
    out = phenotype("CYP2D6", "*4", "*10")
    assert out["activity_score"] == 0.25 and out["phenotype"] == "Intermediate metaboliser"
    assert phenotype("CYP2D6", "*4", "*4")["phenotype"] == "Poor metaboliser"
    assert phenotype("CYP2D6", "*1", "*1x2")["phenotype"] == "Ultrarapid metaboliser"
    assert phenotype("CYP2D6", "1", "41")["activity_score"] == 1.5  # bare numbers accepted


def test_slco1b1_and_tpmt():
    assert phenotype("SLCO1B1", "C", "C")["phenotype"] == "Poor function"
    assert phenotype("tpmt", "*1", "*3A")["phenotype"] == "Intermediate metaboliser"


def test_report_carries_recommendations_and_sources():
    out = report("CYP2C19", "*2", "*2")
    drugs = {r["drug"]: r for r in out["recommendations"]}
    assert "prasugrel" in drugs["Clopidogrel"]["recommendation"].lower()
    assert drugs["Clopidogrel"]["url"].startswith("https://cpicpgx.org/")
    assert len(panel({"CYP2C19": ("*1", "*2"), "SLCO1B1": ("T", "C")})) == 2
    assert {g["gene"] for g in catalogue()} == {"CYP2C19", "CYP2D6", "TPMT", "SLCO1B1"}


def test_unknown_allele_rejected():
    with pytest.raises(BioError):
        phenotype("CYP2C19", "*1", "*99")


# --------------------------------------------------------------------- cli
def test_cli(tmp_path, capsys):
    assert main(["examples", "--out", str(tmp_path)]) == 0
    assert main(["stats", str(tmp_path / "gene.fa")]) == 0
    assert main(["tree", str(tmp_path / "family.fa")]) == 0
    assert main(["variants", str(tmp_path / "reference_cds.fa"), str(tmp_path / "sample_cds.fa")]) == 0
    assert main(["pgx", "CYP2D6", "*4", "*4"]) == 0
    out = capsys.readouterr().out
    assert ";" in out and "frameshift" in out and "Poor metaboliser" in out
    assert main(["pgx", "CYP2D6", "*4", "*99"]) == 2
