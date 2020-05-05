import sys
import argparse

from typing import Mapping
from typing import Sequence, List, Set
from typing import Tuple

from pronto import Ontology, Term

from gffpal.gff import GFF
from gffpal.exceptions import GPInvalidType

"""
Most of this mapping is based on information here:
    https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
    https://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation/
    http://www.insdc.org/insdc-controlled-vocabularies
    https://www.ncbi.nlm.nih.gov/genomes/locustag/Proposal.pdf
    https://www.ncbi.nlm.nih.gov/genbank/evidence/
    http://www.insdc.org/documents/feature-table
    https://www.ncbi.nlm.nih.gov/genbank/collab/db_xref/

Mapping of regulatory elements and repeats is described here:
    https://www.ncbi.nlm.nih.gov/refseq/functionalelements/


Other tools:
    http://genomeannotation.github.io/GAG/
    http://genomeannotation.github.io/annie/
"""


# This is a mapping of SO terms to INDSC terms in descending order of
# specificity. Break after the first subclass_of match.
# Processed is subclass of pseudogene
# Unprocessed is subclass of pseudogene
# Unitary is subclass of unprocessed
# Allelic is subclass of unitary
# Based on: http://www.insdc.org/documents/pseudogene-qualifier-vocabulary
PSEUDOGENE_TYPES = [
    ('SO:0002189', 'allelic'),  # gene
    ('SO:0002200', 'allelic'),  # tRNA
    ('SO:0002196', 'allelic'),  # rRNA
    ('SO:0001759', 'unitary'),  # gene
    ('SO:0002199', 'unitary'),  # tRNA
    ('SO:0002195', 'unitary'),  # rRNA
    ('SO:0001760', 'unprocessed'),  # gene
    ('SO:0002198', 'unprocessed'),  # rRNA
    ('SO:0002194', 'unprocessed'),  # rRNA
    ('SO:0000043', 'processed'),  # gene
    ('SO:0002197', 'processed'),  # tRNA
    ('SO:0002193', 'processed'),  # rRNA
    ('SO:0000336', 'unknown'),  # gene
    ('SO:0000778', 'unknown'),  # tRNA
    ('SO:0000777', 'unknown'),  # rRNA
]


# Repeat info: http://www.insdc.org/controlled-vocabulary-rpttype-qualifier
# Satellite info:
# http://www.insdc.org/controlled-vocabulary-satellite-qualifier

# This is a mapping of SO terms to INDSC terms in descending order of
# specificity. Break after the first subclass_of match.
# Based on: http://www.insdc.org/documents/ncrna-vocabulary
NCRNA_TYPES = [
    ('SO:0001035', 'piRNA'),
    ('SO:0000276', 'miRNA'),
    ('SO:0000646', 'siRNA'),
    ('SO:0002095', 'scaRNA'),
    ('SO:0000013', 'scRNA'),
    ('SO:0000454', 'rasiRNA'),
    ('SO:0000602', 'guide_RNA'),
    ('SO:0000390', 'telomerase_RNA'),
    ('SO:0000385', 'RNase_MRP_RNA'),
    ('SO:0000386', 'RNase_P_RNA'),
    ('SO:0001877', 'lncRNA'),
    ('SO:0000644', 'antisense_RNA'),
    ('SO:0000405', 'Y_RNA'),
    ('SO:0000404', 'vault_RNA'),
    ('SO:0000590', 'SRP_RNA'),
    ('SO:0000274', 'snRNA'),
    ('SO:0000275', 'snoRNA'),
    ('SO:0001998', 'sgRNA'),
    ('SO:0000380', 'hammerhead_ribozyme'),
    ('SO:0000374', 'ribozyme'),
    ('SO:0000588', 'autocatalytically_spliced_intron'),
    ('SO:0001244', 'pre_miRNA'),
    ('SO:0000655', 'other')
]

NCRNA_EXCLUDE_TYPES = [
    "SO:0000253",  # tRNA
    "SO:0000252",  # rRNA
]


TRNA_PRODUCTS = {
    "tRNA": "tRNA-Xxx",
    "alanyl_tRNA": "tRNA-Ala",
    "glutaminyl_tRNA": "tRNA-Gln",
    "prolyl_tRNA": "tRNA-Pro",
    "glutamyl_tRNA": "tRNA-Glu",
    "methionyl_tRNA": "tRNA-Met",
    "asparaginyl_tRNA": "tRNA-Asn",
    "threonyl_tRNA": "tRNA-Thr",
    "glycyl_tRNA": "tRNA-Gly",
    "valyl_tRNA": "tRNA-Val",
    "tyrosyl_tRNA": "tRNA-Tyr",
    "cysteinyl_tRNA": "tRNA-Cys",
    "isoleucyl_tRNA": "tRNA-Ile",
    "seryl_tRNA": "tRNA-Ser",
    "leucyl_tRNA": "tRNA-Leu",
    "selenocysteinyl_tRNA": "tRNA-Sec",
    "tryptophanyl_tRNA": "tRNA-Trp",
    "pyrrolysyl_tRNA": "tRNA-Pyl",
    "lysyl_tRNA": "tRNA-Lys",
    "aspartyl_tRNA": "tRNA-Asp",
    "arginyl_tRNA": "tRNA-Arg",
    "histidyl_tRNA": "tRNA-His",
    "mt_tRNA": "Mitochondrial tRNA",
    "phenylalanyl_tRNA": "tRNA-Phe",
}


def cli_ncbi(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="Input GFF3 file to tidy up. Use '-' for stdin."
    )

    parser.add_argument(
        "--so",
        type=str,
        default="https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/Ontology_Files/so.obo",  # noqa
        help="The sequence ontology OBO file.",
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Output GFF3 file path. Default stdout.",
    )

    parser.add_argument(
        "-l", "--locus-tag",
        type=str,
        default="LCS",
        help="The locus tag prefix to use for genes and locus ids.",
    )

    parser.add_argument(
        "-n", "--number-digits",
        type=int,
        default=5,
        help="",
    )

    parser.add_argument(
        "-p", "--pad-digits",
        type=int,
        default=1,
        help="",
    )

    parser.add_argument(
        "--gene-id-format",
        type=str,
        default="{locus_tag}",
        help="The format for gene feature ids."
    )

    parser.add_argument(
        "--rna-id-format",
        type=str,
        default="{type_lowercase}.{locus_tag}{alpha_index}",
        help="The format for RNA feature ids."
    )

    parser.add_argument(
        "--cds-id-format",
        type=str,
        default="{locus_tag}{parent_alpha_index}",
        help="The format for CDS feature ids."
    )

    parser.add_argument(
        "--other-id-format",
        type=str,
        default="{type_lowercase}{index}.{locus_tag}{parent_alpha_index}",
        help="The format for other feature ids (e.g exons)."
    )

    parser.add_argument(
        "--cds-name-format",
        type=str,
        default="{locus_tag}{parent_alpha_index}"
    )

    parser.add_argument(
        "--transcript-id-format",
        type=str,
        default=(
            "gnl|{locus_tag_prefix}|"
            "{type_lowercase}.{locus_tag}{alpha_index}"
        ),
        help=(
            "What the `transcript_id` field in rna features should be. "
            "The default is based on the NCBIs suggested naming scheme."
        )
    )

    parser.add_argument(
        "--protein-id-format",
        type=str,
        default=(
            "gnl|{locus_tag_prefix}|"
            "{locus_tag}{parent_alpha_index}"
        ),
        help=(
            "What the `protein_id` field in CDS features should be. "
            "The default is based on the NCBIs suggested naming scheme."
        )
    )

    parser.add_argument(
        "--default-rna-product",
        type=str,
        default="hypothetical {type}"
    )

    parser.add_argument(
        "--default-cds-product",
        type=str,
        default="hypothetical protein"
    )

    parser.add_argument(
        "--rna-isoform-product-suffix",
        type=str,
        default="isoform {alpha_index}"
    )

    parser.add_argument(
        "--cds-isoform-product-suffix",
        type=str,
        default="isoform {parent_alpha_index}"
    )

    parser.add_argument(
        "--pseudogene-suffix",
        type=str,
        default="pseudogene"
    )
    return


def add_so_as_ontologies(gff: GFF, name_to_so: Mapping[str, Term]) -> None:
    for record in gff:
        record.add_attributes_if_none()
        assert record.attributes is not None

        if record.type not in name_to_so:
            raise GPInvalidType(record)

        so_term = name_to_so[record.type]
        record.attributes.ontology_term.append(so_term.id)
    return


def add_ncrna_types(
    gff: GFF,
    name_to_so: Mapping[str, Term],
    so: Ontology,
    ncrna_types: Sequence[Tuple[str, str]] = NCRNA_TYPES,
    ncrna_exclude_types: List[str] = NCRNA_EXCLUDE_TYPES,
) -> None:
    ncrna_so_terms = [
        (s, i, set(so[s].subclasses()))
        for s, i
        in ncrna_types
    ]

    ncrna_exclude_terms: Set[str] = set()
    for s in ncrna_exclude_types:
        ncrna_exclude_terms.update(so[s].subclasses())

    for record in gff:
        record.add_attributes_if_none()
        assert record.attributes is not None

        if record.type not in name_to_so:
            raise GPInvalidType(record)

        so_term = name_to_so[record.type]

        if "ncRNA_class" in record.attributes.custom:
            continue
        elif so_term in ncrna_exclude_terms:
            continue

        replacement = None
        for s, i, sub in ncrna_so_terms:
            if so_term in sub:
                record.attributes.custom["ncRNA_class"] = i
                replacement = so[s]
                break

        if replacement is not None:
            if so_term.id not in record.attributes.ontology_term:
                record.attributes.ontology_term.append(so_term.id)

            record.type = replacement.name

    return


def add_pseudogene_types(
    gff: GFF,
    name_to_so: Mapping[str, Term],
    so: Ontology,
    pseudogene_types: Sequence[Tuple[str, str]] = PSEUDOGENE_TYPES,
) -> None:
    pseudogene_so_terms = [
        (s, i, set(so[s].subclasses()))
        for s, i
        in pseudogene_types
    ]

    for record in gff:
        record.add_attributes_if_none()
        assert record.attributes is not None

        if "pseudogene" in record.attributes.custom:
            continue

        if record.type not in name_to_so:
            raise GPInvalidType(record)

        so_term = name_to_so[record.type]
        replacement = None
        for s, i, sub in pseudogene_so_terms:
            if so_term in sub:
                record.attributes.custom["pseudogene"] = i
                record.attributes.custom["pseudo"] = "true"
                replacement = so[s]
                break

        if replacement is not None:
            if so_term.id not in record.attributes.ontology_term:
                record.attributes.ontology_term.append(so_term.id)

            record.type = replacement.name
    return


def replace_genes(gff: GFF, name_to_so: Mapping[str, Term]) -> None:
    gene_subtypes = set(name_to_so["gene"].subclasses())
    gene_subtypes.update(name_to_so["pseudogene"].subclasses())

    gene_subtypes = {
        term.name
        for term
        in name_to_so["gene"].subclasses()
    }

    for record in gff:
        pass
    return


def ncbi(args: argparse.Namespace) -> None:
    gff = GFF.parse(args.infile).break_bubbles()
    so = Ontology.from_obo_library(args.so)

    name_to_so = {
        term.name: term
        for term
        in so.values()
    }

    add_so_as_ontologies(gff, name_to_so)
    add_ncrna_types(gff, name_to_so, so, NCRNA_TYPES)
    add_pseudogene_types(gff, name_to_so, so, PSEUDOGENE_TYPES)

    return
