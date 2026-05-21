#!/usr/bin/env python3

import argparse
from pathlib import Path
from Bio import SeqIO


def validate_fasta(path):
    path = Path(path)
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    ids = [record.id for record in records]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate FASTA headers detected.")

    for record in records:
        if len(record.seq) == 0:
            raise ValueError(f"Empty sequence found for record {record.id}")

    return records


def write_normalized_fasta(records, out_path):
    with open(out_path, "w") as handle:
        for record in records:
            handle.write(f">{record.id}\n")
            handle.write(str(record.seq).strip() + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate and normalize FASTA inputs.")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=False, help="Optional normalized output FASTA")
    args = parser.parse_args()

    records = validate_fasta(args.input)
    print(f"Validated {len(records)} FASTA record(s) from {args.input}.")

    if args.output:
        write_normalized_fasta(records, args.output)
        print(f"Wrote normalized FASTA to {args.output}.")
