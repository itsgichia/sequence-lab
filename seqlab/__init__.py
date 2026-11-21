"""seqlab: small, readable bioinformatics routines, from sequence to tree to prescription.

Pure Python and NumPy. Nothing is fetched from the network.
"""

from .align import align
from .examples import example_family, example_gene, example_variants, examples
from .pgx import GENES, catalogue, panel, phenotype, report
from .sequence import (
    BioError,
    Record,
    clean,
    find_orfs,
    kind_of,
    parse_fasta,
    reverse_complement,
    stats,
    translate,
)
from .tree import jukes_cantor, layout, neighbour_joining, newick, p_distance, tree
from .variants import Variant, call_variants

__version__ = "0.1.0"

__all__ = [
    "BioError",
    "GENES",
    "Record",
    "Variant",
    "align",
    "call_variants",
    "catalogue",
    "clean",
    "example_family",
    "example_gene",
    "example_variants",
    "examples",
    "find_orfs",
    "jukes_cantor",
    "kind_of",
    "layout",
    "neighbour_joining",
    "newick",
    "p_distance",
    "panel",
    "parse_fasta",
    "phenotype",
    "report",
    "reverse_complement",
    "stats",
    "translate",
    "tree",
]
