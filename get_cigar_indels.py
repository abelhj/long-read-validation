import sys
import pysam
from optparse import OptionParser

def print_dels(bamFile, region, minimum):
    if region == "NA":
        iterator = bamFile.fetch()
    else:
        iterator = bamFile.fetch(region=region)
    for read in iterator:
        if not read.is_duplicate and not read.is_secondary:
            contig = read.reference_name
            currentPosition = read.reference_start
            for (operation, length) in read.cigar:
                if operation == 2 and length >= minimum: #deletion
                    print '\t'.join([contig, str(currentPosition-1), str(currentPosition+1), contig, str(currentPosition+length-1), str(currentPosition+length+1), read.query_name, ".", ".", ".", ".", ".", ".", ".", ".", ".", '-'.join(["SVLEN=", str(length)]), "DEL"])
                if operation == 1 and length >= minimum: #insertion
                    print '\t'.join([contig, str(currentPosition-1), str(currentPosition+1), contig, str(currentPosition-1), str(currentPosition+1), read.query_name, ".", ".", ".", ".", ".", ".", ".", ".", ".", '='.join(["SVLEN", str(length)]), "INS"])
                if operation !=1 and operation != 4 and operation !=5 and operation !=6: #insertions and clips don't advance the position
                    currentPosition = currentPosition + length

def main():
    usage = "%prog -a <bam file> -r <region> -m <minimumLength>\nVersion: 0.1\nAuthor: Allison Regier\n"
    parser = OptionParser(usage)
    parser.add_option("-a", dest="bamFile", metavar="FILE", help="BAM file")
    parser.add_option("-r", dest="region", metavar="String", help="region", default="NA")
    parser.add_option("-m", dest="minimum", metavar="Int", help="minimum deletion length to print", default=50)

    (opts, args) = parser.parse_args()

    if opts.bamFile is None:
        parser.print_help()
        print
    else:
        bamFile = pysam.AlignmentFile(opts.bamFile, "rb")
        print_dels(bamFile, opts.region, opts.minimum)

if __name__ == "__main__":
    sys.exit(main())
