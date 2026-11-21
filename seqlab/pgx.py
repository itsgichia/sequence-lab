"""Pharmacogenomics: from a genotype to a metaboliser phenotype to a dosing recommendation.

This is the personalised-medicine end of the toolkit. It encodes, in a small
and inspectable form, the logic that clinical pharmacogenomics services apply:

1. Each star allele (or SNP genotype) of a pharmacogene has a *function*
   (normal, decreased, no function, increased) and, for CYP2D6, an *activity
   value*.
2. A diplotype (two alleles) maps to a *phenotype*: poor, intermediate, normal,
   rapid or ultrarapid metaboliser.
3. Each phenotype carries a *recommendation* for selected drugs, condensed from
   the published CPIC guidelines (Clinical Pharmacogenetics Implementation
   Consortium, https://cpicpgx.org/guidelines/).

The tables are deliberately small and the text is paraphrased. They exist to
show how the pipeline works, not to guide treatment. Allele function
assignments can change as evidence accumulates and real reports consider
far more alleles, copy-number variants and phasing. Nothing here is medical
advice.
"""

from __future__ import annotations

from .sequence import BioError

PM, IM, NM, RM, UM = (
    "Poor metaboliser",
    "Intermediate metaboliser",
    "Normal metaboliser",
    "Rapid metaboliser",
    "Ultrarapid metaboliser",
)
NORMAL_FN, DECREASED_FN, NO_FN, INCREASED_FN, POOR_FN = "normal", "decreased", "no function", "increased", "poor"

GENES: dict[str, dict] = {
    "CYP2C19": {
        "name": "Cytochrome P450 2C19",
        "role": "Activates the antiplatelet prodrug clopidogrel; clears proton pump inhibitors and several antidepressants.",
        "input": "star",
        "alleles": {
            "*1": NORMAL_FN,
            "*2": NO_FN,
            "*3": NO_FN,
            "*17": INCREASED_FN,
        },
        "notes": "*2 (c.681G>A, rs4244285) and *3 (c.636G>A, rs4986893) create splicing and nonsense changes that abolish activity. *17 (rs12248560) is a promoter variant that increases expression.",
        "drugs": {
            "Clopidogrel": {
                "guideline": "CPIC guideline for CYP2C19 and clopidogrel (2022 update)",
                "url": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/",
                PM: "Avoid standard-dose clopidogrel. Poor metabolisers form little of the active metabolite and have higher rates of stent thrombosis and cardiovascular events. Use prasugrel or ticagrelor at standard dose if there is no contraindication.",
                IM: "Avoid standard-dose clopidogrel where an alternative is appropriate. Reduced active-metabolite formation; use prasugrel or ticagrelor at standard dose if there is no contraindication.",
                NM: "Standard dose clopidogrel.",
                RM: "Standard dose clopidogrel.",
                UM: "Standard dose clopidogrel.",
            },
            "Omeprazole": {
                "guideline": "CPIC guideline for CYP2C19 and proton pump inhibitors (2020)",
                "url": "https://cpicpgx.org/guidelines/cpic-guideline-for-proton-pump-inhibitors-and-cyp2c19/",
                PM: "Standard starting dose for acute treatment. For chronic therapy beyond 12 weeks, consider a 50% dose reduction once efficacy is established, because exposure is much higher.",
                IM: "Standard starting dose. For chronic therapy, consider a 50% reduction once efficacy is established.",
                NM: "Standard starting dose. Consider increasing by 50 to 100% for Helicobacter pylori eradication or erosive oesophagitis.",
                RM: "Increase the starting dose by 50 to 100% for H. pylori eradication or erosive oesophagitis; monitor for efficacy.",
                UM: "Increase the starting dose by 50 to 100% for H. pylori eradication or erosive oesophagitis; monitor for efficacy.",
            },
        },
    },
    "CYP2D6": {
        "name": "Cytochrome P450 2D6",
        "role": "Metabolises about a quarter of commonly used drugs, including codeine (activation to morphine), tramadol, ondansetron, tamoxifen and many antidepressants.",
        "input": "star",
        "activity": {
            "*1": 1.0, "*2": 1.0, "*3": 0.0, "*4": 0.0, "*5": 0.0, "*6": 0.0,
            "*9": 0.5, "*10": 0.25, "*17": 0.5, "*29": 0.5, "*41": 0.5,
            "*1x2": 2.0, "*2x2": 2.0, "*4x2": 0.0,
        },
        "alleles": {
            "*1": NORMAL_FN, "*2": NORMAL_FN, "*3": NO_FN, "*4": NO_FN, "*5": NO_FN, "*6": NO_FN,
            "*9": DECREASED_FN, "*10": DECREASED_FN, "*17": DECREASED_FN, "*29": DECREASED_FN, "*41": DECREASED_FN,
            "*1x2": INCREASED_FN, "*2x2": INCREASED_FN, "*4x2": NO_FN,
        },
        "notes": "Phenotype comes from the activity score, the sum of the two allele values: 0 is a poor metaboliser, 0.25 to 1 intermediate, 1.25 to 2.25 normal, above 2.25 ultrarapid (2019 CPIC and DPWG consensus). *4 is the common no-function allele in Europeans, *10 (decreased) is common in East Asians and *17 and *29 in people of African ancestry. xN marks a gene duplication.",
        "drugs": {
            "Codeine": {
                "guideline": "CPIC guideline for CYP2D6, OPRM1 and COMT and opioid therapy (2021)",
                "url": "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/",
                PM: "Avoid codeine: too little morphine is formed for pain relief. Use a non-opioid or an opioid that does not depend on CYP2D6 activation (not tramadol).",
                IM: "Use the label dose. If there is no response, consider a non-codeine opioid; avoid tramadol as the alternative.",
                NM: "Use the label dose.",
                UM: "Avoid codeine: rapid conversion to morphine risks toxicity, including respiratory depression. Choose a non-opioid or an opioid that is not a CYP2D6 substrate.",
            },
            "Ondansetron": {
                "guideline": "CPIC guideline for CYP2D6 and ondansetron and tropisetron (2017)",
                "url": "https://cpicpgx.org/guidelines/guideline-for-ondansetron-and-tropisetron-and-cyp2d6-genotype/",
                PM: "Standard dose. Limited evidence of altered response.",
                IM: "Standard dose.",
                NM: "Standard dose.",
                UM: "Select an alternative antiemetic not predominantly metabolised by CYP2D6, such as granisetron, because of reduced efficacy.",
            },
        },
    },
    "TPMT": {
        "name": "Thiopurine S-methyltransferase",
        "role": "Inactivates thiopurines (azathioprine, mercaptopurine, thioguanine) used in leukaemia, inflammatory bowel disease and transplantation.",
        "input": "star",
        "alleles": {
            "*1": NORMAL_FN,
            "*2": NO_FN,
            "*3A": NO_FN,
            "*3B": NO_FN,
            "*3C": NO_FN,
            "*4": NO_FN,
        },
        "notes": "*3A (two missense changes, Ala80Pro and Tyr240Cys) is the most common no-function allele in Europeans; *3C (Tyr240Cys alone) is more common in African and East Asian populations. About 1 in 300 people carry two no-function alleles.",
        "drugs": {
            "Azathioprine": {
                "guideline": "CPIC guideline for TPMT and NUDT15 and thiopurine dosing (2018 update)",
                "url": "https://cpicpgx.org/guidelines/guideline-for-thiopurines-and-tpmt/",
                PM: "Drastically reduce the dose: for malignancy, reduce the daily dose tenfold and give three times weekly instead of daily; for non-malignant conditions consider an alternative non-thiopurine immunosuppressant. Risk of life-threatening myelosuppression at standard doses.",
                IM: "Start at 30 to 80% of the normal dose and titrate on blood counts after 2 to 4 weeks.",
                NM: "Start at the normal dose and adjust on blood counts after 2 weeks.",
            },
        },
    },
    "SLCO1B1": {
        "name": "Organic anion transporting polypeptide 1B1",
        "role": "Carries statins from blood into liver cells. Reduced transport raises plasma statin levels and the risk of muscle toxicity.",
        "input": "snp",
        "variant": {"rsid": "rs4149056", "cdna": "c.521T>C", "protein": "p.Val174Ala", "ref": "T", "alt": "C"},
        "alleles": {"T": NORMAL_FN, "C": DECREASED_FN},
        "notes": "A single missense variant explains most of the clinically used variation; c.521C is also the defining change of the *5 and *15 alleles. Entered as a genotype of two bases rather than star alleles.",
        "drugs": {
            "Simvastatin": {
                "guideline": "CPIC guideline for SLCO1B1, ABCG2, CYP2C9 and statins (2022)",
                "url": "https://cpicpgx.org/guidelines/cpic-guideline-for-statins/",
                "Poor function": "Avoid simvastatin. Prescribe an alternative statin at a dose with a lower myopathy risk (for example rosuvastatin or pravastatin), and consider the indication and other risk factors.",
                "Decreased function": "Prescribe an alternative statin, or if simvastatin is used keep the dose at 20 mg a day or less and monitor for muscle symptoms.",
                "Normal function": "Prescribe the desired starting dose and adjust on efficacy and tolerability.",
            },
        },
    },
}

_CPIC_FUNCTION_PHENOTYPE = {
    (NO_FN, NO_FN): PM,
    (NO_FN, DECREASED_FN): PM,
    (DECREASED_FN, DECREASED_FN): IM,
    (NO_FN, NORMAL_FN): IM,
    (DECREASED_FN, NORMAL_FN): IM,
    (NO_FN, INCREASED_FN): IM,
    (DECREASED_FN, INCREASED_FN): IM,
    (NORMAL_FN, NORMAL_FN): NM,
    (NORMAL_FN, INCREASED_FN): RM,
    (INCREASED_FN, INCREASED_FN): UM,
}


def _normalise_allele(gene: dict, allele: str) -> str:
    a = allele.strip().upper()
    if gene["input"] == "star":
        if not a.startswith("*"):
            a = "*" + a
        a = a.replace("XN", "x2").replace("X2", "x2")
    if a not in gene["alleles"]:
        known = ", ".join(gene["alleles"])
        raise BioError(f"Unknown allele '{allele}'. This table knows: {known}.")
    return a


def phenotype(gene_symbol: str, allele_1: str, allele_2: str) -> dict:
    """Map a diplotype to a metaboliser (or transporter function) phenotype."""
    gene = GENES.get(gene_symbol.upper())
    if gene is None:
        raise BioError(f"Unknown gene '{gene_symbol}'. Available: {', '.join(GENES)}.")
    a1, a2 = _normalise_allele(gene, allele_1), _normalise_allele(gene, allele_2)
    f1, f2 = gene["alleles"][a1], gene["alleles"][a2]
    out = {
        "gene": gene_symbol.upper(),
        "diplotype": f"{a1}/{a2}",
        "allele_functions": {a1: f1, a2: f2} if a1 != a2 else {a1: f1},
        "activity_score": None,
    }
    if "activity" in gene:
        score = gene["activity"][a1] + gene["activity"][a2]
        out["activity_score"] = score
        if score == 0:
            pheno = PM
        elif score <= 1.0:
            pheno = IM
        elif score <= 2.25:
            pheno = NM
        else:
            pheno = UM
    elif gene["input"] == "snp":
        n_dec = (f1 == DECREASED_FN) + (f2 == DECREASED_FN)
        pheno = ("Normal function", "Decreased function", "Poor function")[n_dec]
    else:
        key = tuple(sorted((f1, f2), key=[NO_FN, DECREASED_FN, NORMAL_FN, INCREASED_FN].index))
        pheno = _CPIC_FUNCTION_PHENOTYPE[key]
    out["phenotype"] = pheno
    return out


def report(gene_symbol: str, allele_1: str, allele_2: str) -> dict:
    """Phenotype plus the drug recommendations for it."""
    out = phenotype(gene_symbol, allele_1, allele_2)
    gene = GENES[out["gene"]]
    recs = []
    for drug, table in gene["drugs"].items():
        text = table.get(out["phenotype"])
        if text is None:  # phenotype absent from the table (e.g. no RM for CYP2D6)
            text = table.get(NM, "No specific recommendation in this table.")
        recs.append({"drug": drug, "recommendation": text, "guideline": table["guideline"], "url": table["url"]})
    out["recommendations"] = recs
    out["gene_name"] = gene["name"]
    out["gene_role"] = gene["role"]
    out["notes"] = gene["notes"]
    return out


def panel(genotypes: dict[str, tuple[str, str]]) -> list[dict]:
    """Reports for several genes at once, e.g. {"CYP2C19": ("*1", "*2"), "SLCO1B1": ("T", "C")}."""
    return [report(g, a1, a2) for g, (a1, a2) in genotypes.items()]


def catalogue() -> list[dict]:
    """The gene table in a form a user interface can render (alleles, drugs, notes)."""
    out = []
    for symbol, g in GENES.items():
        out.append(
            {
                "gene": symbol,
                "name": g["name"],
                "role": g["role"],
                "input": g["input"],
                "variant": g.get("variant"),
                "alleles": [{"allele": a, "function": f, "activity": g.get("activity", {}).get(a)} for a, f in g["alleles"].items()],
                "drugs": [{"drug": d, "guideline": t["guideline"], "url": t["url"]} for d, t in g["drugs"].items()],
                "notes": g["notes"],
            }
        )
    return out


__all__ = ["GENES", "phenotype", "report", "panel", "catalogue"]
