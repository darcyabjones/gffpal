# GFFPal

A quick and flexible toolset/python library for working with GFF files.

I'm often needing to convert formats into or from GFF3.

Genometools and Aegean are useful tools for some aspects of GFF manipulation (especially tidying and conversion to GTF).
But they often aren't great at handling non-specifcation compliant GFF3s (i.e. CDS attached directly to genes, "match" style GFFs, or augustus hints GFFs).
Also I wanted something that I could use to select by ID (instead of tricky `awk` or `grep -f` commands).

This is intended for my own use and I'm motivated by ease of use rather than computational effeciency (hence Python).
But if you have thoughts, contributions, or whatever, please raise an issue to get in touch.


## Available scripts

All scripts are subcommands under the gffpal command.
To view a list of subcommands type `gffpal -h`.

#### `gffpal hints`

Converts gffs in conventional gene formats into "hints" for the [augustus](https://github.com/Gaius-Augustus/Augustus) gene prediction toolset.
This is similar in spirit to the `align2hints.pl` script in [BRAKER](https://github.com/Gaius-Augustus/BRAKER/tree/master/scripts) but we make no effort to add new features (e.g. start, stop), or parse the various non-standard formats.
It's a simple conversion of the type columns, grouping genes using the `group` attribute, and adding priorities etc.
To add additional features (e.g. start, stop, UTR etc), I recommend `canon-gff3` from [AEGeAn](https://aegean.readthedocs.io/en/stable/).

Basic usage:

```
gffpal hints \
  --outfile "my_hints.gff3" \
  --source "M" \
  --priority 4 \
  --group-level "mRNA" \
  --cds "CDSpart" \
  --cds-trim 12 \
  genes_to_create_hints_from.gff3
```


#### `gffpal expandcds`

A utility to fix wonky CDS features in GFF3 files.
The GFF3 standard says that start and stop codons should be included in the CDS
ranges but many tools exclude the stop codon.

This script simply expands the CDS by 3 bases to include the start and/or codons.
Optionally if you provide the GFFs corresponding genome fasta we will print out warnings for all of the non-standard start/stop codons according to some translation table.
So you can evaluate if you're actually doing the right thing :)

NOTE: at the moment the start/stop extension and extraction is really stupid.
We simply bump the number from the relevant CDS section by three.
This could give incorrect results for cases where there is an intron within the start/stop codon.
Future versions might use exon information to handle this, but for now it's an acceptable limitation for my use case (Fixing spaln output prior to augustus hint generation).

Basic usage:

```
gffpal expandcds \
  --outfile "expanded.gff3" \
  --infasta "my_genome.fasta" \
  --gencode 1 \
  --warnings "non_standard.gff3" \
  "my_genome.gff3"
```


#### `gffpal rnammer2gff`

Converts [RNAmmer](http://www.cbs.dtu.dk/services/RNAmmer/) gff2 output to gff3 with the correct "type"
column set according to the [Sequence Ontology](http://www.sequenceontology.org/).

Basic usage:

```
rnammer -S euk -m lsu,ssu,tsu -gff "rnammer.gff2" my.fasta
gffpal rnammer2gff -o rnammer.gff3 -k euk rnammer.gff2
```


#### `gffpal trnascan2gff`

Converts [tRNAscan](http://lowelab.ucsc.edu/tRNAscan-SE/) output to gff3 with the correct "types"
column set according to the [Sequence Ontology](http://www.sequenceontology.org/) and including codon/secondary structure information in the attributes..

Basic usage:

```
tRNAscan-SE -E -o "trnascan.txt" -f "trnascan_ss.txt" "my.fasta"

gffpal trnascan2gff -o "trnascan.gff3" trnascan.txt trnascan_ss.txt
```
