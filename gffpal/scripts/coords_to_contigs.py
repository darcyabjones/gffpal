import sys
import argparse

from math import log10, ceil

from typing import Sequence, List
from typing import Iterable, Iterator
from typing import Set
from typing import Tuple
from typing import Dict, Mapping
from typing import Optional

from copy import copy

from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import align

from intervaltree import IntervalTree, Interval

from gffpal.gff import GFF3Record, Strand, Phase
from gffpal.gff import GFF3Attributes
from gffpal.parsers.coords import Coords


def cli_coords(parser):
    parser.add_argument(
        "coords",
        type=argparse.FileType('r'),
        help="The mummer coords alignment file.",
    )

    parser.add_argument(
        "scaffolds",
        type=argparse.FileType('r'),
        help="The scaffolds fasta file.",
    )

    parser.add_argument(
        "contigs",
        type=argparse.FileType('r'),
        help="The contigs fasta file.",
    )

    parser.add_argument(
        "-e", "--extend",
        action="store_true",
        default=False,
        help=(
            "Try to extend the contigs to include adjacent non-N sequence. "
            "Trim regions of Ns from the ends of sequences, and create new "
            "contigs for short regions not covered by N-stretches or contigs."
        )
    )

    parser.add_argument(
        "-p", "--min-pid",
        type=float,
        default=0.95,
        help="The minimum percent identity [0-1]",
    )

    parser.add_argument(
        "-q", "--min-qcov",
        type=float,
        default=0.95,
        help="The minimum contig coverage to be considered a match [0-1]",
    )

    parser.add_argument(
        "-a", "--min-aln",
        type=float,
        default=125,
        help=(
            "The minimum alignment length to consider as a match in bp. "
            "This should be something just at or below your minimum contig "
            "size."
        ),
    )

    parser.add_argument(
        "--mum-source",
        type=str,
        default="MUMmer",
        help="What should be in the source column for the coord matches."
    )

    parser.add_argument(
        "--source",
        type=str,
        default="gffpal",
        help="What should be in the source column for the predicted contigs."
    )

    parser.add_argument(
        "-n", "--nstretch",
        type=int,
        default=10,
        help=(
            "We check stretches of Ns to see if contigs overlap them and "
            "if the whole genome is covered by contigs."
        )
    )

    parser.add_argument(
        "-r", "--no-print-nstretch",
        dest="print_nstretch",
        action="store_false",
        default=True,
        help=(
            "Don't print the locations of long stretches of N in the contig "
            "gff. These features will all have the type 'gap' (SO:0000730); "
            "however, in many cases the length of the gap will be unknown. "
            "In the absence of a SO term for gap with unknown length, i'm "
            "using this. Please resolve in your own way if you need to."
        )
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Where to write the tidied gff to."
    )

    parser.add_argument(
        "-m", "--matches",
        action="store_true",
        default=False,
        help="Write the un-processed mummer matches out to the gff3 file too."
    )
    return


def find_n_stretches(
    scaffolds: Mapping[str, SeqRecord],
    n_stretch: int,
) -> Dict[str, IntervalTree]:
    itrees = dict()

    for scaffold, seq in scaffolds.items():
        itrees[scaffold] = IntervalTree()

        n_count = 0
        for i, base in enumerate(seq):
            if base.upper() == "N":
                n_count += 1

            elif n_count >= n_stretch:
                interval = Interval(i - n_count, i)
                itrees[scaffold].add(interval)
                n_count = 0

            else:
                n_count = 0

        if n_count >= n_stretch:
            interval = Interval(i - n_count, i)
            itrees[scaffold].add(interval)

    return itrees


def find_identical_seqs(
    contigs: Mapping[str, SeqRecord],
    scaffolds: Mapping[str, SeqRecord]
) -> Set[Tuple[str, str, Strand]]:
    """ Sometimes mummer misses identical sequences because
    there are no unique matches, (the MUMs). This stupid match method
    based on checksums just catches those extra bits.
    """
    from Bio.SeqUtils.CheckSum import seguid

    same: Dict[str, List[Tuple[str, Strand]]] = defaultdict(list)
    for contig, contig_seq in contigs.items():
        checksum = seguid(contig_seq.seq)
        rchecksum = seguid(contig_seq.seq.reverse_complement())

        same[checksum].append((contig, Strand.PLUS))
        same[rchecksum].append((contig, Strand.MINUS))

    matches = set()
    for scaffold, scaffold_seq in scaffolds.items():
        checksum = seguid(scaffold_seq.seq)
        contig_matches = same[checksum]

        if len(contig_matches) == 0:
            continue

        # Just in case the checksum is a collision.
        filtered_matches = []
        for contig, strand in contig_matches:
            contig_seq = contigs[contig].seq
            if strand == Strand.MINUS:
                contig_seq = contig_seq.reverse_complement()

            if str(contig_seq) == str(scaffold_seq.seq):
                filtered_matches.append((contig, strand))

        matches.add((scaffold, contig, strand))

    return matches


def add_identical_seqs_to_coords(
    matches: Set[Tuple[str, str, Strand]],
    coords: Iterable[Coords],
    contigs: Mapping[str, SeqRecord],
    scaffolds: Mapping[str, SeqRecord],
) -> List[Coords]:
    out: List[Coords] = []

    for scaffold, contig, strand in matches:
        scaffold_len = len(scaffolds[scaffold])
        contig_len = len(contigs[contig])
        out.append(Coords(
            0, scaffold_len,
            0, contig_len,
            strand,
            scaffold_len, contig_len,
            100.0,
            scaffold_len, contig_len,
            100.0, 100.0,
            scaffold, contig,
        ))

    for coord in coords:
        tup = (coord.ref, coord.query, coord.strand)

        if tup not in matches:
            out.append(coord)
    return out


def filter_coords(
    coords: Iterable[Coords],
    min_pid: float = 0.95,
    min_qcov: float = 0.6,
    min_aln: int = 100,
) -> Iterator[Coords]:
    for coord in coords:
        if coord.qalnlen < min_aln:
            continue

        if (coord.pid / 100) < min_pid:
            continue

        if (coord.qcov / 100) < min_qcov:
            continue

        yield coord
    return


def coords_to_itrees(coords: Iterable[Coords]) -> Dict[str, IntervalTree]:
    itrees: Dict[str, IntervalTree] = defaultdict(IntervalTree)
    for coord in coords:
        itrees[coord.ref].add(coord.as_interval())
    return itrees


def intersection(left, right):
    lstart = min([left.begin, left.end])
    lend = max([left.begin, left.end])
    rstart = min([right.begin, right.end])
    rend = max([right.begin, right.end])

    start = max([lstart, rstart])
    end = min([lend, rend])
    return Interval(start, end)


def find_len_non_overlap(interval: Interval, itree: IntervalTree) -> int:
    overlaps = IntervalTree(itree.overlap(interval))
    overlaps.merge_overlaps()

    len_overlap = sum([
        intersection(interval, o).length()
        for o
        in overlaps
    ])

    return interval.length() - len_overlap


def find_best_match_in_cluster(
    interval: Interval,
    itree: IntervalTree,
    coords_mapping: Mapping[Coords, Interval],
    min_non_overlap: int
) -> Tuple[Optional[Interval], List[Interval]]:
    cant_match = []
    matches = []

    for coord in interval.data:
        coord_interval = coords_mapping[coord]
        len_non_overlap = find_len_non_overlap(
            coord_interval,
            itree
        )
        if len_non_overlap < min_non_overlap:
            cant_match.append(coord_interval)
        else:
            matches.append((len_non_overlap, coord_interval))

    if len(matches) == 0:
        return None, cant_match

    matches.sort(
        key=lambda t: (t[0], t[1].data.pid, t[1].data.qcov),
        reverse=True
    )

    return matches[0][1], cant_match


def select_best_alignments(
    itrees: Mapping[str, IntervalTree],
    coords: Iterable[Coords],
    min_non_overlap: int = 80,
) -> None:
    """ MUTATES ITREES """

    done = set()
    for scaffold, itree in itrees.items():
        for interval in itree:
            done.add(interval.data)

    coords = list(coords)
    coords_itrees = coords_to_itrees(c for c in coords if c not in done)

    for scaffold, coords_itree in coords_itrees.items():
        while len(coords_itree) > 0:
            coords_itree_copy = copy(coords_itree)
            coords_mapping = {i.data: i for i in coords_itree_copy}
            coords_itree_copy.merge_overlaps(
                data_reducer=lambda x, y: x + [y],
                data_initializer=[]
            )

            for interval in coords_itree_copy:
                best, cant_match = find_best_match_in_cluster(
                    interval,
                    itrees[scaffold],
                    coords_mapping,
                    min_non_overlap,
                )

                if best is not None:
                    itrees[scaffold].add(best)

                for coord_interval in cant_match:
                    coords_itree.remove(coord_interval)
    return


def filter_contained_contigs(itrees: Mapping[str, IntervalTree]) -> None:
    """ Remove any contigs that are completely contained within other contigs.
    """
    for scaffold, itree in itrees.items():
        to_drop = set()
        for interval in itree:
            for contained in itree.envelop(interval):
                if contained != interval:
                    to_drop.add(contained)
        for interval in to_drop:
            itree.remove(interval)
    return


def filter_nstretches(
    itrees: Mapping[str, IntervalTree],
    nstretches: Mapping[str, IntervalTree],
    min_non_overlap: int,
) -> None:
    """
    Remove contigs without much going on outside N stretches.
    """

    for scaffold, itree in nstretches.items():

        # Loop through all of the potential breaks
        for nstretch in itree:
            # Find "contigs" that overlap the potential break
            # We do this in sorted order, from smallest to largest alignment
            # that means shorter ones are removed first
            contigs = sorted(
                itrees[scaffold].overlap(nstretch),
                key=lambda x: x.length()
            )

            to_drop = set()

            # Loop through the contigs to test.
            for contig in contigs:
                # Find if they overlap with any other n stretches.
                n_overlaps = nstretches[scaffold].overlap(contig)

                # Get an intervaltree of all contigs overlapping this one.
                contig_overlaps = IntervalTree(
                    itrees[scaffold].overlap(contig)
                )

                # Remove all of the n-chunks from the intervals.
                # Note the "Coords" is still duplicated in the data attribute
                for n_overlap in n_overlaps:
                    contig_overlaps.chop(n_overlap.begin, n_overlap.end)

                # Get the intervals that aren't the overlap under
                # consideration.
                contig_overlaps_itree = IntervalTree(
                    o
                    for o
                    in contig_overlaps
                    if o.data != contig.data
                )
                contig_overlaps_itree.merge_overlaps()

                # Get the fragments of the overlap under consideration
                contig_itree = IntervalTree(
                    o
                    for o
                    in contig_overlaps
                    if o.data == contig.data
                )

                # For each of the fragments, find how many new Non-N bases it
                # contributes to the contigging.
                len_non_overlap = sum([
                    find_len_non_overlap(f, contig_overlaps_itree)
                    for f
                    in contig_itree
                ])

                # Remove the contig if it doesn't cut the muster
                if len_non_overlap < min_non_overlap:
                    to_drop.add(contig)

            for contig in to_drop:
                itrees[scaffold].remove(contig)
    return


def filter_non_informative(
    itrees: Mapping[str, IntervalTree],
    min_non_overlap: int,
) -> None:
    """
    This is to handle the case of a non-informative hit overlapping
    two other contigs but without overlapping an nstretch
    eg.

    <-----------> <---------->
              !--------!
    """

    for scaffold, itree in itrees.items():
        itree_copy = copy(itree)
        coords_mapping = {i.data: i for i in itree_copy}
        itree_copy.merge_overlaps(
            data_reducer=lambda x, y: x + [y],
            data_initializer=[]
        )
        for cluster in itree_copy:
            cluster_itree = IntervalTree(
                coords_mapping[c]
                for c
                in cluster.data
            )
            while len(cluster_itree) > 0:
                matches = []
                for contig in cluster_itree:
                    contig_itree = IntervalTree(
                        c
                        for c
                        in cluster_itree
                        if c != contig
                    )

                    len_non_overlap = find_len_non_overlap(
                        contig,
                        contig_itree
                    )

                    matches.append((len_non_overlap, contig))

                matches.sort(
                    key=lambda t: (t[0], t[1].data.pid, t[1].data.qcov)
                )
                len_non_overlap, worst = matches[0]
                if len_non_overlap < min_non_overlap:
                    cluster_itree.remove(worst)
                    itrees[scaffold].remove(worst)
                else:
                    # By definition we can't remove any more
                    break
    return


def decide_if_complete_match(left: Coords, right: Coords) -> List[Coords]:

    left_complete = left.qend == left.qlen
    right_complete = right.qstart == 0

    if left_complete and right_complete:
        return [left, right]
    elif left_complete:
        return [left]
    elif right_complete:
        return [right]
    else:
        return [left, right]


def find_best_nstretch(
    interval: Interval,
    left: Coords,
    right: Coords,
    nstretches: Sequence[Interval],
) -> Interval:

    best = set()
    best_dist = 100000000
    for nstretch in nstretches:
        left_dist_to_qend = left.qlen - left.qend
        right_dist_to_qstart = right.qstart

        left_dist = min([
            abs(nstretch.begin - left.rstart) + left_dist_to_qend,
            abs(nstretch.end - left.rstart) + left_dist_to_qend,
        ])

        if left_dist < best_dist:
            best = {nstretch}
            best_dist = left_dist
        elif left_dist == best_dist:
            best.add(nstretch)

        right_dist = min([
            abs(nstretch.begin - right.rend) + right_dist_to_qstart,
            abs(nstretch.end - right.rend) + right_dist_to_qstart,
        ])

        if right_dist < best_dist:
            best = {nstretch}
            best_dist = right_dist
        elif right_dist == best_dist:
            best.add(nstretch)

    if len(best) > 1:
        max_length = max(b.length() for b in best)
        best = {b for b in best if b.length() == max_length}

    assert len(best) == 1
    return best.pop()


def extract_scaffold_seq(seq: SeqRecord, interval: Interval) -> str:
    return str(seq.seq[interval.begin: interval.end])


def extract_contig_seqs(
    contigs: Mapping[str, SeqRecord],
    interval: Interval,
    pad: int = 5,
) -> List[Tuple[Coords, str]]:
    out = []

    for coord in interval.data:
        # Rstart etc is the offset within the alignment.
        # We pad the sequences a bit because we don't have gap info
        rstart = max([0, interval.begin - coord.rstart - pad])
        rend = min([coord.rlen, (interval.end - coord.rstart) + pad])
        seqrecord = contigs[coord.query]

        if coord.strand == Strand.MINUS:
            seq = seqrecord.seq[coord.qstart: coord.qend].reverse_complement()
        else:
            seq = seqrecord.seq[coord.qstart: coord.qend]

        out.append((coord, str(seq)[rstart: rend]))

    return out


def align_intersections(
    contigs: List[Tuple[Coords, str]],
    scaffold_seq: str,
) -> Coords:
    max_score = 0
    best = None

    for coord, contig_seq in contigs:
        # This scoring is based on what minimap uses
        score = align.localms(
            scaffold_seq,
            contig_seq,
            2,  # match
            -4,  # mismatch
            -5,  # gapopen
            -2,  # gapextend
            penalize_end_gaps=False,
            score_only=True
        )
        if score > max_score:
            best = coord
            max_score = score

    assert best is not None
    return best


def split_at_nstretch(
    interval: Interval,
    left: Coords,
    right: Coords,
    nstretches: Sequence[Interval],
) -> List[Interval]:
    best_nstretch = find_best_nstretch(interval, left, right, nstretches)
    itree = IntervalTree([interval])
    itree.chop(best_nstretch.begin, best_nstretch.end)
    split_intervals = list(itree)

    # This should be 0 (if the whole intersection is N), 1 (if the nstretch
    # buts one of the ends), or 2 (if internal split).
    assert len(split_intervals) <= 2, split_intervals
    split_intervals_with_data = []
    for i in split_intervals:
        assert not ((i.begin == interval.begin) and (i.end == interval.end))

        data = None
        if i.begin == interval.begin:
            data = left
        elif i.end == interval.end:
            data = right

        assert data is not None
        split_intervals_with_data.append(Interval(
            i.begin,
            i.end,
            data
        ))
    return split_intervals_with_data


def split_overlaps(
    itrees: Mapping[str, IntervalTree],
    nstretches: Mapping[str, IntervalTree],
    scaffolds: Mapping[str, SeqRecord],
    contigs: Mapping[str, SeqRecord],
) -> None:
    for scaffold, itree in itrees.items():
        itree.split_overlaps()
        itree.merge_overlaps(
            data_reducer=lambda x, y: x + [y],
            data_initializer=[]
        )

        for interval in itree:
            if len(interval.data) == 1:
                continue

            # By this point we should have filtered out anything that could
            # give more than 2 alignments at a position.
            assert len(interval.data) == 2, interval

            data = sorted(interval.data, key=lambda x: x.rstart)
            left = data[0]
            right = data[1]

            # This should be true if it's sorted.
            assert left.rend == interval.end, left
            assert right.rstart == interval.begin, right

            n_overlaps = nstretches[scaffold].overlap(interval)
            if len(n_overlaps) > 0:
                split = split_at_nstretch(interval, left, right, n_overlaps)
                itree.remove(interval)
                for i in split:
                    itree.add(i)

                continue

            interval_data = decide_if_complete_match(left, right)
            if len(interval_data) == 1:
                itree.remove(interval)
                itree.add(Interval(
                    interval.begin,
                    interval.end,
                    interval_data
                ))
                continue

            scaffold_seq = extract_scaffold_seq(scaffolds[scaffold], interval)
            contig_seqs = extract_contig_seqs(contigs, interval)

            best = align_intersections(contig_seqs, scaffold_seq)
            itree.remove(interval)
            itree.add(Interval(
                interval.begin,
                interval.end,
                [best]
            ))

    return


def merge_adjacent_intervals(itree: IntervalTree) -> List[Interval]:
    intervals = sorted(itree, key=lambda x: (x.begin, x.end))
    out = list()

    previous = None
    for interval in intervals:
        if previous is None:
            previous = interval
        elif ((interval.begin == previous.end) and
              (interval.data == previous.data)):
            # Extend the interval!
            previous = Interval(previous.begin, interval.end, previous.data)
        elif interval.begin >= previous.end:
            out.append(previous)
            previous = interval
        else:
            raise AssertionError(
                "Somehow we got two overlapping intervals at this stage."
                f"left: {previous}, right: {interval}."
            )

    if previous is not None:
        out.append(previous)

    return out


def extend_contigs(
    itrees: Mapping[str, IntervalTree],
    nstretches: Mapping[str, IntervalTree],
    scaffolds: Mapping[str, SeqRecord]
) -> None:
    for scaffold, seq in scaffolds.items():
        new_contigs: List[Interval] = list()
        contigs = sorted(itrees[scaffold], key=lambda x: x.begin)
        nstretches_itree = nstretches[scaffold]

        for i in range(len(contigs)):

            if i == 0:
                last_contig_end = 0
            else:
                last_contig_end = new_contigs[-1].end

            if i == len(contigs) - 1:
                next_contig_start = len(seq)
            else:
                next_contig_start = contigs[i + 1].begin

            contig = contigs[i]

            if contig.begin > last_contig_end:
                ns = nstretches_itree.overlap(last_contig_end, contig.begin)
                if len(ns) > 0:
                    n = sorted(ns, key=lambda x: x.end)[-1]
                    # If an stretch is actually at the junction, this could
                    # bump the contig start back.
                    contig = Interval(n.end, contig.end, contig.data)
                elif last_contig_end == 0:
                    contig = Interval(0, contig.end, contig.data)
                else:
                    contig = Interval(last_contig_end, contig.end, contig.data)

            if contig.end < next_contig_start:
                ns = nstretches_itree.overlap(contig.end, next_contig_start)
                if len(ns) > 0:
                    n = sorted(ns, key=lambda x: x.begin)[0]
                    # If an stretch is actually at the junction, this could
                    # bump the contig start back.
                    contig = Interval(contig.begin, n.begin, contig.data)
                elif next_contig_start == len(seq):
                    contig = Interval(contig.begin,
                                      next_contig_start, contig.data)
                else:
                    halfway = (contig.end +
                               ceil((next_contig_start - contig.end) / 2))
                    contig = Interval(contig.begin, halfway, contig.data)

            new_contigs.append(contig)

        itrees[scaffold] = IntervalTree(new_contigs)
    return


def trim_contigs(
    itrees: Mapping[str, IntervalTree],
    scaffolds: Mapping[str, SeqRecord],
) -> None:
    for scaffold, seq in scaffolds.items():
        new_contigs = list()

        contigs = sorted(itrees[scaffold], key=lambda x: x.begin)
        last_end = 0

        for i in range(len(contigs)):
            contig = contigs[i]
            assert contig.begin >= last_end

            last_end = contig.end

            j = contig.begin
            while (j < len(seq)) and (seq.seq[j] == "N"):
                j += 1

            if j != contig.begin:
                contig = Interval(j, contig.end, contig.data)

            k = contig.end - 1
            while (k >= 0) and (seq.seq[k] == "N"):
                k -= 1

            if k != contig.end - 1:
                contig = Interval(contig.begin, k + 1, contig.data)

            if len(contig) > 1:
                new_contigs.append(contig)
        itrees[scaffold] = IntervalTree(new_contigs)
    return


def find_remaining(
    itrees: Mapping[str, IntervalTree],
    nstretches: Mapping[str, IntervalTree],
    scaffolds: Mapping[str, SeqRecord],
) -> None:
    for scaffold, seq in scaffolds.items():
        contigs = itrees[scaffold]
        nstretch = nstretches[scaffold]

        # This is just to remove the data from the intervals.
        # Having data prevents them from being removed with difference.
        intervals = [Interval(i.begin, i.end) for i in contigs]
        intervals.extend(Interval(i.begin, i.end) for i in nstretch)
        covered = IntervalTree(intervals)
        # Strict=false means that adjacent but non-overlapping
        # will also be merged.
        covered.merge_overlaps(strict=False)

        remaining = IntervalTree([Interval(0, len(seq))]) | covered
        remaining.split_overlaps()
        remaining.difference_update(covered)

        itrees[scaffold].update(remaining)
    return


def coord_to_interval(
    interval: Interval,
    index: int,
    ndigits: int,
    source: str,
) -> GFF3Record:
    assert len(interval.data) == 1, interval
    gffrecord = interval.data[0].as_gffrecord()
    gffrecord.start = interval.begin
    gffrecord.end = interval.end
    gffrecord.type = "contig"
    gffrecord.attributes.id = f"contig{index:0>{ndigits}}"

    # Because we split some of them, i want the alignment info in the
    # attributes.
    gffrecord.attributes.custom.update({
        "scaffold_alignment_start": interval.data[0].rstart,
        "scaffold_alignment_end": interval.data[0].rend,
    })
    gffrecord.attributes.note.append(
        "Contigs aligned to scaffolds using MUMmer"
    )

    gffrecord.source = source
    return gffrecord


def remaining_to_gffrecord(
    seqid: str,
    interval: Interval,
    index: int,
    ndigits: int,
    source: str,
) -> GFF3Record:
    return GFF3Record(
        seqid,
        source,
        "contig",
        interval.begin,
        interval.end,
        None,
        Strand.UNSTRANDED,
        Phase.NOT_CDS,
        GFF3Attributes(id=f"contig{index:0>{ndigits}}")
    )


def gap_to_gffrecord(seqid: str, interval: Interval) -> GFF3Record:
    return GFF3Record(
        seqid,
        "gffpal",
        "gap",
        interval.begin,
        interval.end,
        None,
        Strand.UNSTRANDED,
        Phase.NOT_CDS,
        GFF3Attributes()
    )


def coords(args: argparse.Namespace) -> None:  # noqa
    contig_seqs = SeqIO.to_dict(SeqIO.parse(args.contigs, "fasta"))
    scaffold_seqs = SeqIO.to_dict(SeqIO.parse(args.scaffolds, "fasta"))

    identical_seqs = find_identical_seqs(contig_seqs, scaffold_seqs)

    coords = add_identical_seqs_to_coords(
        identical_seqs,
        Coords.from_file(args.coords),
        contig_seqs,
        scaffold_seqs,
    )

    gffrecords = list()
    if args.matches:
        for coord in filter_coords(
            coords,
            args.min_pid,
            args.min_qcov,
            args.min_aln
        ):
            gffrecords.append(coord.as_gffrecord(source=args.mum_source))

    itrees: Dict[str, IntervalTree] = defaultdict(IntervalTree)

    # We assign contig positions in 3 stages with lower stringency.
    # This is just so that the best matches take priority.

    # SELECT_BEST_ALIGNMENTS MUTATES ITREES
    select_best_alignments(
        itrees,
        filter_coords(coords, 0.99, 0.99, args.min_aln),
        args.min_aln,
    )

    select_best_alignments(
        itrees,
        filter_coords(coords, 0.95, args.min_qcov, args.min_aln),
        args.min_aln,
    )

    filter_contained_contigs(itrees)
    nstretches = find_n_stretches(scaffold_seqs, args.nstretch)
    filter_nstretches(itrees, nstretches, args.min_aln)
    filter_non_informative(itrees, args.min_aln)
    split_overlaps(itrees, nstretches, scaffold_seqs, contig_seqs)

    # Here we flatten the split contigs, merging butted ones.
    total_n_contigs = 0
    contigs_intervals = dict()
    for scaffold, seq in sorted(
        scaffold_seqs.items(),
        key=lambda x: len(x[1]),
        reverse=True,
    ):
        contigs_intervals[scaffold] = merge_adjacent_intervals(
            itrees.get(scaffold, IntervalTree())
        )
        total_n_contigs += len(contigs_intervals[scaffold])

    if args.extend:
        extend_contigs(contigs_intervals, nstretches, scaffold_seqs)
        find_remaining(contigs_intervals, nstretches, scaffold_seqs)
        trim_contigs(contigs_intervals, scaffold_seqs)

    # This is just for padding the contig id numbers.
    ndigits = ceil(log10(total_n_contigs))

    index = 0
    for scaffold, contigs in sorted(
        contigs_intervals.items(),
        key=lambda t: t[0]
    ):
        for contig in sorted(contigs, key=lambda x: x.begin):
            if contig.data is None:
                gffrecord = remaining_to_gffrecord(scaffold, contig,
                                                   index, ndigits, args.source)
            else:
                gffrecord = coord_to_interval(contig, index,
                                              ndigits, args.source)

            gffrecords.append(gffrecord)
            index += 1

    if args.print_nstretch:
        for scaffold, nstretch in nstretches.items():
            for n in nstretch:
                gffrecords.append(gap_to_gffrecord(scaffold, n))

    print("##gff-version 3", file=args.outfile)
    for scaffold, sseq in sorted(
        scaffold_seqs.items(),
        key=lambda x: len(x[1]),
        reverse=True,
    ):
        print(
            f"##sequence-region   {scaffold} 1 {len(sseq)}",
            file=args.outfile
        )

    for gffrecord in sorted(
        gffrecords,
        key=lambda x: (x.seqid, x.start, x.end, x.type)
    ):
        print(gffrecord, file=args.outfile)
    return
