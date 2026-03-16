#!/usr/bin/env python3
"""Edit a SAM file to make it graph-friendly.

samtools will find reads that overlap a region,
but I need reads that are contained within.
Also I need to shift the coordinates.

So this script will edit the start coordinate to
be graph-friendly, and drop any reads which are
not completely contained within the region.
"""

import argparse # Command line argument parsing
import re # CIGAR-parsing regex

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')

# Operations that consume reference
REF_CONSUME = {'M', 'D', 'N', '=', 'X'}

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Filter and shift coordinates in a headerless SAM file.'
    )
    parser.add_argument('sam', help='Headerless SAM file')
    parser.add_argument('-s', '--start', type=int, help='Start coordinate')
    parser.add_argument('-e', '--end', type=int, help='End coordinate')
    return parser.parse_args()

def cigar_ref_length(cigar: str) -> int:
    """Return reference length consumed by a CIGAR string."""
    length = 0
    for n, op in CIGAR_RE.findall(cigar):
        if op in REF_CONSUME:
            length += int(n)
    return length

if __name__ == '__main__':
    args = parse_args()

    with open(args.sam) as alignments:
        for line in alignments:
            line = line.strip()
            if not line:
                # Empty line; probably at the end of the file
                continue

            # https://samtools.github.io/hts-specs/SAMv1.pdf
            fields = line.split('\t')
            pos = int(fields[3])
            cigar = fields[5]

            ref_len = cigar_ref_length(cigar)
            end_pos = pos + ref_len - 1

            # Drop if alignment extends past end
            if end_pos > args.end:
                continue

            # +1 because the dummy node's N at the start
            new_pos = pos - args.start + 1

            # Drop if new start would be invalid
            if new_pos < 0:
                continue

            fields[3] = str(new_pos)
            print('\t'.join(fields))