#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: bash 04_generate_cuttag_normalized_bigwigs.sh <manifest.tsv> <output_directory> [threads]" >&2
  exit 1
fi

manifest=$1
output_dir=$2
threads=4
if [[ $# -eq 3 ]]; then
  threads=$3
fi

if [[ ! -f "$manifest" ]]; then
  echo "Manifest not found: $manifest" >&2
  exit 1
fi
if ! command -v bamCoverage >/dev/null 2>&1; then
  echo "bamCoverage was not found. Install deepTools 3.5.1 before running this script." >&2
  exit 1
fi
if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
  echo "threads must be a positive integer" >&2
  exit 1
fi

mkdir -p "$output_dir"

header=$(head -n 1 "$manifest")
expected_header=$'sample_id\tcondition\tbam_path\toutput_filename'
if [[ "$header" != "$expected_header" ]]; then
  echo "Expected tab-separated header: $expected_header" >&2
  exit 1
fi

while IFS=$'\t' read -r sample_id condition bam_path output_filename; do
  [[ -z "$sample_id$condition$bam_path$output_filename" ]] && continue
  [[ "$sample_id" == \#* ]] && continue

  if [[ -z "$sample_id" || -z "$condition" || -z "$bam_path" || -z "$output_filename" ]]; then
    echo "Incomplete manifest row for sample: $sample_id" >&2
    exit 1
  fi
  if [[ ! -f "$bam_path" ]]; then
    echo "BAM not found for $sample_id: $bam_path" >&2
    exit 1
  fi
  if [[ "$output_filename" == */* || "$output_filename" != *.bw ]]; then
    echo "output_filename must be a .bw basename: $output_filename" >&2
    exit 1
  fi

  output_path="$output_dir/$output_filename"
  echo "Generating $output_path from $sample_id ($condition)" >&2
  bamCoverage \
    --bam "$bam_path" \
    --outFileName "$output_path" \
    --outFileFormat bigwig \
    --normalizeUsing CPM \
    --ignoreForNormalization chrX \
    --binSize 50 \
    --numberOfProcessors "$threads"
done < <(tail -n +2 "$manifest")

echo "CPM-normalized CUT&Tag BigWig files written to $output_dir" >&2
