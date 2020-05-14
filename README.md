# GFFPal

A quick and flexible toolset/python library for working with GFF files.

I'm often needing to convert formats into or from GFF3.

Genometools and Aegean are useful tools for some aspects of GFF manipulation (especially tidying and conversion to GTF).
But they often aren't great at handling non-specifcation compliant GFF3s (i.e. CDS attached directly to genes, "match" style GFFs, or augustus hints GFFs).
Also I wanted something that I could use to select by ID (instead of tricky `awk` or `grep -f` commands).

This is intended for my own use and I'm motivated by ease of use rather than computational effeciency (hence Python).
But if you have thoughts, contributions, or whatever, please raise an issue to get in touch.


## Installing

At the moment things tend to just be written directly to the master github branch and I don't do much versioning.

Until that changes, i'd suggest installing from here :)
If you want versioned point releases let me know.

If you have `git` and `python3` installed you can do:

```
python3 -m pip install git+https://github.com/darcyabjones/gffpal.git
```

If you aren't using a virtualenv or conda env, i'd suggest using the `--user` option to avoid messing with your root packages.

```
python3 -m pip install --user git+https://github.com/darcyabjones/gffpal.git
```

There are ways to ask it to install specific commits, or to use the master zip if you don't have `git` installed.
Create an issue if you feel it would help you.


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


#### `gffpal exonerate2gff`

Converts exonerates weird gff2 format to a GFF3 with unique ids for each block.
This is a lossy process, some of the ids that exonerate generates are lost (but the query protein is retained).
It's only known to work on the protein to genome alignments, so others might have issues.

The source is always set to "exonerate".
"exon" features are used as "CDSs" since protein alignments should only give that and they have the match statistics.

Basic usage:

```
gffpal exonerate2gff -o "exonerate.gff3" exonerate.gff2
```


#### `gffpal coord2contig`


Coord to contig attempts to find the best non-overlapping set of contigs given
a nucmer coords alignment to the scaffolds of the same assembly.

Essentially it just iteratively selects alignments based on their length, % identity, contig coverage, and how much they overlap with other contig alignments.
The remaining overlapping alignments are split at the intersection of the two based on the following:

1. If there is a stretch of N-s in the intersection, split the alignment intersection at the one closest to one of the alignment ends.
2. If one (and only one) of the contig alignments goes to the end of the contig, then the whole alignment intersection gets assigned to that one and the two contigs will "butt" against each other without an N-stretch.
3. If there is an n-stretch immediately next to the intersection, the one that is unbroken gets the region.
4. Otherwise, extract the sequences corresponding to the intersection for each contig, locally align it to the scaffold section, and take the highest scoring match.
   In the case of a tie, the left-most contig will win.

Optionally, with the `--extend` parameter, the program will try to extend the "contigs" to fill space between contigs or the ends of sequences, trim any contig ends at N's, and add contigs to fill in any stretches of unassigned genome that don't overlap an N-stretch longer than specified by the `--nstretch` parameter. Effectively this becomes similar to the old strategy of splitting the genome at N-stretches, but you won't necessarily split contigs that have good alignments crossing them and you'll catch cases of butting contigs with no Ns separating them.

Note that the SPAdes assembler doesn't give a nice 1-1 contig tiling path, some of the contigs overlap, and sometimes it will add some extra sequence between the contigs.
The extend option can be useful for finding something approximating a tiling path.

To use, align your scaffolds and contigs using Nucmer and get the coords file.
I used MUMmer version 4.0.0beta2.


```
nucmer --breaklen 400 --threads 4 \
  --delta genome.delta \
  scaffolds.fasta \
  contigs.fasta

show-coords -T -c -l -o genome.delta > genome.coords
```

Then you can run coord2contig like so:

```
gffpal coord2contig -o contig.gff3 genome.coords scaffolds.fasta contigs.fasta
```


You can also get the good matches from the coords file as a gff (before trying to resolve overlaps) using the `--matches` flag.

Mummer alignments will be given the type `nucleotide_match` (SO:0000347).
Contigs will be given the `contig` type (SO:0000149), and n-stretches will have the `gap` type (SO:0000730).


Use `gffpal coord2contig --help` for full options.
