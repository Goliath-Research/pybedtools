_bedtools_path = ""
_R_path = ""

tempfile_prefix = "pybedtools."
tempfile_suffix = ".tmp"

# Checking for BEDTools will happen when creating the first BedTool; other
# checks happen at first use (BAM object creation; tabix-ing a BedTool)
_bedtools_installed = False
_R_installed = False
_v_2_15_plus = False
_v_2_27_plus = False
bedtools_version = []

KEEP_TEMPFILES = False
_DEBUG = True

# Ceilings (inclusive, in bytes) for gzip-compressed interval files that
# pybedtools reads itself (BedTool iteration, head, saveas, cat, to_dataframe,
# and IntervalFile). Going past one raises pybedtoolsError. Set a ceiling to
# None to disable it. Files handed to BEDTools programs are read by bedtools.
gzip_max_compressed_bytes = 128 * 1024
gzip_max_decompressed_bytes = 256 * 1024
gzip_max_line_bytes = 64 * 1024

# Check calls against these names to only allow calls to known BEDTools
# programs (basic security)
#
_prog_names = {
    # Genome arithmetic
    "intersectBed": "intersect",
    "windowBed": "window",
    "closestBed": "closest",
    "coverageBed": "coverage",
    "mapBed": "map",
    "genomeCoverageBed": "genomecov",
    "mergeBed": "merge",
    "clusterBed": "cluster",
    "complementBed": "complement",
    "subtractBed": "subtract",
    "slopBed": "slop",
    "shiftBed": "shift",
    "flankBed": "flank",
    "sortBed": "sort",
    "randomBed": "random",
    "shuffleBed": "shuffle",
    "annotateBed": "annotate",
    "spacing": "spacing",
    # multi-way
    "multiIntersectBed": "multiinter",
    "unionBedGraphs": "unionbedg",
    # PE
    "pairToBed": "pairtobed",
    "pairToPair": "pairtopair",
    # format conversion
    "bamToBed": "bamtobed",
    "bedToBam": "bedtobam",
    "bedpeToBam": "bedpetobam",
    "bed12ToBed6": "bed12tobed6",
    "bamToFastq": "bamtofastq",
    # fasta
    "fastaFromBed": "getfasta",
    "maskFastaFromBed": "maskfasta",
    "nucBed": "nuc",
    # bam-centric
    "multiBamCov": "multicov",
    "tagBam": "tag",
    # stats
    "jaccard": "jaccard",
    "reldist": "reldist",
    # misc
    "getOverlap": "overlap",
    "bedToIgv": "igv",
    "linksBed": "links",
    "windowMaker": "makewindows",
    "groupBy": "groupby",
    "expandCols": "expand",
    "sample": "sample",
    "fisher": "fisher",
    "split": "split",
}

_old_names = list(_prog_names.keys())
_new_names = list(_prog_names.values())

_column_names = {
    "bed": [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRgb",
        "blockCount",
        "blockSizes",
        "blockStarts",
    ],
    "gff": [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attributes",
    ],
}
