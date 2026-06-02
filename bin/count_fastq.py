#!/usr/bin/python3

import argparse
import gzip
import logging
import sys
import csv

logging.basicConfig(level=logging.INFO, format="%(levelname)s : %(message)s")


class ReadReportParser(argparse.ArgumentParser):
    def error(self, msg):
        self.print_help()
        sys.stderr.write(f"\nERROR: {msg}\n")
        sys.exit(1)

if __name__ == "__main__":
    parser = ReadReportParser(
        description="Generate a CSV report with read counts"
    )

    parser.add_argument(
        "--sample_id",
        help="Sample ID for the report"
    )
    parser.add_argument(
        "--raw_reads",
        nargs="+",
        required=True,
        help="Raw read files"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV file"
    )

    args = parser.parse_args()

FIELDS = [
    "id",
    "count1", "count2"
    ]

def write_csv(path, row):
    """Write a single-row CSV."""
    try:
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerow(row)
    except Exception as e:
        logging.error(f"Cannot write CSV {path}: {e}")
        sys.exit(1)

    logging.info(f"Report written: {path}")

def count_reads(path):
    """Count reads in FASTQ."""
    opener = gzip.open if path.endswith(".gz") else open

    try:
        with opener(path, "rt") as f:
            lines = sum(1 for _ in f)
    except Exception as e:
        logging.error(f"Cannot read FASTQ {path}: {e}")
        sys.exit(1)

    if lines % 4:
        logging.error(f"FASTQ not divisible by 4 (possible malformed): {path}")

    return lines // 4

def process(raw, output, sample_id):
    """Main processing function. """

    logging.info("Counting reads")
    raw_counts = [count_reads(f) for f in raw]

    row = {
        "id":              sample_id,
        "count1":             raw_counts[0],
        "count2":             raw_counts[1]
    }

    write_csv(output, row)

process(args.raw_reads, args.output, args.sample_id)
