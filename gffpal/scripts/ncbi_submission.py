import sys
import argparse

from pronto import Ontology

from gffpal import GFF


pseudogene_types = [
    {
        "type": "processed",
        "description": "The pseudogene has arisen by reverse transcription of a mRNA into cDNA, followed by reintegration into the genome. Therefore, it has lost any intron/exon structure, and it might have a pseudo-polyA-tail.",
        "so_term": "SO:0000043",
        "priority": 2
    },  # subclass of "pseudogene"
    {
        "type": "unprocessed",
        "description": "The pseudogene has arisen from a copy of the parent gene by duplication followed by accumulation of random mutation. The changes, compared to their functional homolog, include insertions, deletions, premature stop codons, frameshifts and a higher proportion of non-synonymous versus synonymous substitutions.",
        "so_term": "SO:0001760",
        "priority": 2
    },  # subclass of "pseudogene"
    {
        "type": "unitary",
        "description": "The pseudogene has no parent. It is the original gene, which is functional in some species but disrupted in some way (indels, mutation, recombination) in another species or strain.",
        "so_term": "SO:0001759",
        "priority": 3

    },  # Unitary is subclass of unprocessed
    {
        "type": "allelic",
        "description": "A (unitary) pseudogene that is stable in the population but importantly it has a functional alternative allele also in the population. i.e., one strain may have the gene, another strain may have the pseudogene. MHC haplotypes have allelic pseudogenes.",
        "so_term": "SO:0002189",
        "priority": 4
    },  # Allelic is subclass of unitary.
    {
        "type": "unknown",
        "description": "The submitter does not know the method of pseudogenisation.",
        "so_term": "SO:0000336",
        "priority": 1
    }
]


ncrna_types = [
    {
        "type": "antisense_RNA",
        "description": "RNA molecule that is transcribed from the opposite strand of DNA compared to another transcript.",
        "so_term": "SO:0000644",
        "priority": 2,
    },
    {
        "type": "autocatalytically_spliced_intron",
        "description": "self-splicing intron.",
        "so_term": "SO:0000588",
        "priority": 1,
    },
    {
        "type": "ribozyme",
        "description": "ribonucleic acid enzyme, RNA molcule that can catalyse specific biochemical reactions.",
        "so_term": "SO:0000374",
        "priority": 1,
    },
    {
        "type": "hammerhead_ribozyme",
        "description": "Small catalytic RNA motif that catalyzes a self-cleavage reaction. Its name comes from its secondary structure which resembles a carpenter's hammer. The hammerhead ribozyme is involved in the replication of some viroid and satellite RNAs.",
        "so_term": "SO:0000380",
        "priority": 1,
    },
    {
        "type": "lncRNA",
        "description": "Long non-coding RNA; such molecules are generally defined as having a length greater than 200bp and do not fit into any other ncRNA class.",
        "so_term": "SO:0001877",
        "priority": 2,
    },
    {
        "type": "RNase_P_RNA",
        "description": "RNA component of Ribonuclease P (RNase P), a ubiquitous endoribonuclease.",
        "so_term": "SO:0000386",
        "priority": 2,
    },
    {
        "type": "RNase_MRP_RNA",
        "description": "RNA molecule essential for the catalytic activity of RNase MRP, an enzymatically active ribonucleoprotein with two distinct roles in eukaryotes. In mitochondria it plays a direct role in the initiation of mitochondrial DNA replication, while in the nucleus it is involved in precursor rRNA processing.",
        "so_term": "SO:0000385",
        "priority": 2,
    },
    {
        "type": "telomerase_RNA",
        "description": "RNA component of telomerase, a reverse transcriptase that synthesises telomeric DNA.",
        "so_term": "SO:0000390",
        "priority": 2,
    },
    {
        "type": "guide_RNA",
        "description": "Short 3'-uridylated RNA that can form a duplex with a stretch of mature edited mRNA.",
        "so_term": "SO:0000602",
        "priority": 2,
    },
    {
        "type": "sgRNA",
        "description": "A small RNA oligo, typically about 20 bases, that guides the cas nuclease to a target DNA sequence in the CRISPR/cas mutagenesis method.",
        "so_term": "SO:0001998",
        "priority": 1,
    },
    {
        "type": "rasiRNA",
        "description": "Small interfering RNA of length between 17 and 28 nucleotides, derived from the transcript of a repetitive element.",
        "so_term": "SO:0000454",
        "priority": 2,
    },
    {
        "type": "scRNA",
        "description": "Small cytoplasmic RNA; any one of several small cytoplasmic RNA molecules present in the cytoplasm and (sometimes) nucleus of a eukaryote.",
        "so_term": "SO:0000013",
        "priority": 2,
    },
    {
        "type": "scaRNA",
        "description": "An ncRNA, specific to the Cajal body, that has been demonstrated to function as a guide RNA in the site-specific synthesis of 2'-O-ribose-methylated nucleotides and pseudouridines in the RNA polymerase II-transcribed U1, U2, U4 and U5 spliceosomal small nuclear RNAs (snRNAs).",
        "so_term": "SO:0002095",
        "priority": 2,
    },
    {
        "type": "siRNA",
        "description": "SO:0000646",
        "so_term": "Small RNA molecule that is the product of a longer exogenous or endogenous double stranded RNA, which is either a bimolecular duplex or very long hairpin, processed (via the Dicer pathway) such that numerous siRNAs accumulate from both strands of the double stranded RNA. sRNAs trigger the cleavage of their target molecules.",
        "priority": 3,
    },
    {
        "type": "pre_miRNA",
        "description": "The 60-70 nucleotide region remaining after Drosha processing of a microRNA primary transcript, where this region folds back upon itself to form a hairpin structure from which a mature microRNA is processed.",
        "so_term": "SO:0001244",
        "priority": 1,
    },
    {
        "type": "miRNA",
        "description": "Small, ~22-nt, RNA molecule, termed microRNA, produced from precursor molecules that can form local hairpin structures, which ordinarily are processed (via the Dicer pathway) such that a single miRNA molecule accumulates from one arm of a hairpin precursor molecule. MicroRNAs may trigger the cleavage of their target molecules or act as translational repressors.",
        "so_term": "SO:0000276",
        "priority": 3,
    },
    {
        "type": "piRNA",
        "description": "Small RNA molecule, termed Piwi-interacting RNA, expressed in testes and forming RNA-protein complex with Piwi protein; purification of these complexes has revealed that Piwi-interacting RNA oligonucleotides are approximately 24-32 nucleotides long",
        "so_term": "SO:0001035",
        "priority": 3,
    },
    {
        "type": "snoRNA",
        "description": "Small [nucleolar] RNA molecules involved in modification and processing of ribosomal RNA or transfer RNA; found in archaea and in eukaryotic species where they are often localized to the nucleolus, but are not necessarily nucleolar-specific, e.g., some subsets may be scaRNAs that are localized to the Cajal body.",
        "so_term": "SO:0000275",
        "priority": 2,
    },
    {
        "type": "snRNA",
        "description": "Small nuclear RNA molecules involved in pre-mRNA splicing and processing.",
        "so_term": "SO:0000274",
        "priority": 2,
    },
    {
        "type": "SRP_RNA",
        "description": "Signal recognition particle, a universally conserved ribonucleoprotein involved in the co-translational targeting of proteins to membranes.",
        "so_term": "SO:0000590",
        "priority": 2,
    },
    {
        "type": "vault_RNA",
        "description": "RNA component of the vault ribonuceoprotein, a complex which consists of a major vault protein (MVP), two minor vault proteins (VPARP and TEP1), and several small RNA molecules and has been suggested to be involved in drug resistance.",
        "so_term": "SO:0000404",
        "priority": 2,
    },
    {
        "type": "Y_RNA",
        "description": "Component of the Ro ribonucleoprotein particle (Ro RNP), in association with Ro60 and La proteins.",
        "so_term": "SO:0000405",
        "priority": 2,
    },
    {
        "type": "other",
        "description": "SO:0000655",
        "so_term": "ncRNA_class not included in any other term",
        "priority": 1,
    },
]


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
        "--mrna-id-format",
        type=str,
        default="{parent_gene}_mRNA{index}",
        help="The format for mRNA feature ids."
    )

    parser.add_argument(
        "--cds-id-format",
        type=str,
        default="{parent_mrna}_CDS",
        help="The format for mRNA feature ids."
    )

    parser.add_argument(
        "--exon-id-format",
        type=str,
        default="{parent_mrna}_exon{index}",
        help="The format for mRNA feature ids."
    )

    parser.add_argument(
        "--cds-name-format",
        type=str,
        default="{parent_gene}_protein{parent_mrna_index}"
    )

    parser.add_argument(
        "--transcript-id-format",
        type=str,
        default="gnl|{locus_tag_prefix}|gene{locus_tag_index}_{type}{index}",
    )

    parser.add_argument(
        "--protein-id-format",
        type=str,
        default="gnl|{locus_tag_prefix}|gene{locus_tag_index}_protein{parent_mrna_index}",
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
    return


def ncbi(args: argparse.Namespace) -> None:
    gff = GFF.parse(args.infile)
    for line in gff
    return
