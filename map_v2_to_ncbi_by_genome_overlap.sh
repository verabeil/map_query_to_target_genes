#!/bin/bash
#SBATCH -o %x.out
#SBATCH -e %x.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=80G
#SBATCH --time=6:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vbeilins@caltech.edu
#SBATCH --job-name=map_v2_ncbi_loci
#SBATCH -A carnegie_poc

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:-10}"

CONDA_INIT="/home/vbeilins/miniconda3/bin/activate"
CONDA_ENV="v2_ncbi_map"

QUERY_LABEL="v2"
QUERY_GENOME="Lachesis_assembly.fasta.gz"
QUERY_GTF="eupsc_models_v2.2.gtf"
QUERY_TRANSCRIPTS=""
QUERY_GTF_EXPLICIT="0"

TARGET_LABEL="ncbi"
TARGET_GENOME="GCF_053919585.1_ASM5391958v1_genomic.fna.gz"
TARGET_GTF="ncbi.gtf"

OUTDIR=""

# Transcript-to-genome alignment filters.
MIN_MAPQ="20"
MIN_QCOV="0.50"

# Annotation overlap filters for assigning a query transcript to a target gene.
MIN_EXON_OVERLAP_BASES="50"
MIN_EXON_OVERLAP_FRAC="0.10"
MIN_BEST_MARGIN="0.10"

usage() {
  cat <<'EOF'
Usage:
  map_v2_to_ncbi_by_genome_overlap.sh [options]

Maps genes from a query annotation to genes in a target annotation by:
  1. extracting query transcripts from the query genome,
  2. aligning them to the target genome,
  3. assigning each query transcript to the target gene whose exons overlap most.

If --query-transcripts is provided, the script skips transcript extraction and uses
that transcript FASTA directly. This is useful for transcriptome assemblies such
as Trinity output.

Options:
  --query-label NAME       Label for query assembly/annotation. Default: v2
  --query-genome FILE      Query genome FASTA, optionally .gz.
  --query-gtf FILE         Query annotation GTF/GFF.
  --query-transcripts FILE Query transcript FASTA, optionally .gz. Skips gffread.
  --target-label NAME      Label for target assembly/annotation. Default: ncbi
  --target-genome FILE     Target genome FASTA, optionally .gz.
  --target-gtf FILE        Target annotation GTF/GFF.
  --outdir DIR             Output directory. Default: QUERY_to_TARGET_genome_overlap
  --threads N              Threads. Default: SLURM_CPUS_PER_TASK or 10
  --min-mapq N             Minimum transcript-genome MAPQ. Default: 20
  --min-qcov FLOAT         Minimum query transcript coverage. Default: 0.50
  --min-exon-bases N       Minimum exon-overlap bases. Default: 50
  --min-exon-frac FLOAT    Minimum fraction of aligned query bases overlapping target exons. Default: 0.10
  --min-best-margin FLOAT  Minimum margin between best and second-best target gene. Default: 0.10
  --conda-init FILE        Conda activation script. Default: /home/vbeilins/miniconda3/bin/activate
  --conda-env NAME         Conda env to activate. Default: v2_ncbi_map
  --no-conda               Do not activate conda inside the script.
  -h, --help               Show this help.

Default inputs are set for the current V2-to-NCBI run.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --query-label) QUERY_LABEL="$2"; shift 2 ;;
    --query-genome) QUERY_GENOME="$2"; shift 2 ;;
    --query-gtf) QUERY_GTF="$2"; QUERY_GTF_EXPLICIT="1"; shift 2 ;;
    --query-transcripts) QUERY_TRANSCRIPTS="$2"; shift 2 ;;
    --target-label) TARGET_LABEL="$2"; shift 2 ;;
    --target-genome) TARGET_GENOME="$2"; shift 2 ;;
    --target-gtf) TARGET_GTF="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
    --min-qcov) MIN_QCOV="$2"; shift 2 ;;
    --min-exon-bases) MIN_EXON_OVERLAP_BASES="$2"; shift 2 ;;
    --min-exon-frac) MIN_EXON_OVERLAP_FRAC="$2"; shift 2 ;;
    --min-best-margin) MIN_BEST_MARGIN="$2"; shift 2 ;;
    --conda-init) CONDA_INIT="$2"; shift 2 ;;
    --conda-env) CONDA_ENV="$2"; shift 2 ;;
    --no-conda) CONDA_ENV=""; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$OUTDIR" ]]; then
  OUTDIR="${QUERY_LABEL}_to_${TARGET_LABEL}_genome_overlap"
fi

if [[ -n "$QUERY_TRANSCRIPTS" ]]; then
  required_inputs=("$QUERY_TRANSCRIPTS" "$TARGET_GENOME" "$TARGET_GTF")
  if [[ "$QUERY_GTF_EXPLICIT" == "1" ]]; then
    required_inputs+=("$QUERY_GTF")
  fi
elif [[ -n "$QUERY_GENOME" && -n "$QUERY_GTF" ]]; then
  required_inputs=("$QUERY_GENOME" "$QUERY_GTF" "$TARGET_GENOME" "$TARGET_GTF")
else
  echo "Provide either --query-transcripts, or both --query-genome and --query-gtf." >&2
  exit 1
fi

for input in "${required_inputs[@]}"; do
  if [[ ! -s "$input" ]]; then
    echo "Missing or empty input file: $input" >&2
    exit 1
  fi
done

if [[ -n "$CONDA_ENV" ]]; then
  source "$CONDA_INIT"
  conda activate "$CONDA_ENV"
fi

if [[ -z "$QUERY_TRANSCRIPTS" ]]; then
  QUERY_TRANSCRIPTS="$OUTDIR/${QUERY_LABEL}.transcripts.fa"
fi
SAM="$OUTDIR/${QUERY_LABEL}_transcripts_to_${TARGET_LABEL}_genome.sam"

if [[ -n "$QUERY_TRANSCRIPTS" && "$QUERY_GTF_EXPLICIT" == "0" && "$QUERY_TRANSCRIPTS" != "$OUTDIR/${QUERY_LABEL}.transcripts.fa" ]]; then
  USE_QUERY_GTF="0"
else
  USE_QUERY_GTF="1"
fi

export QUERY_LABEL TARGET_LABEL QUERY_GTF TARGET_GTF OUTDIR SAM USE_QUERY_GTF
export MIN_MAPQ MIN_QCOV MIN_EXON_OVERLAP_BASES MIN_EXON_OVERLAP_FRAC MIN_BEST_MARGIN

mkdir -p "$OUTDIR"

if [[ -n "$QUERY_GENOME" && -n "$QUERY_GTF" && "$QUERY_TRANSCRIPTS" == "$OUTDIR/${QUERY_LABEL}.transcripts.fa" ]]; then
  echo "Extracting ${QUERY_LABEL} transcript FASTA..."
  gffread "$QUERY_GTF" -g "$QUERY_GENOME" -w "$QUERY_TRANSCRIPTS"
else
  echo "Using provided ${QUERY_LABEL} transcript FASTA: $QUERY_TRANSCRIPTS"
fi

echo "Aligning ${QUERY_LABEL} transcripts to the ${TARGET_LABEL} genome..."
minimap2 -ax splice -t "$THREADS" --secondary=no \
  "$TARGET_GENOME" \
  "$QUERY_TRANSCRIPTS" \
  > "$SAM"

echo "Assigning ${QUERY_LABEL} transcripts to overlapping ${TARGET_LABEL} genes..."
python3 - <<'PY'
import bisect
import os
import re
from collections import defaultdict

outdir = os.environ["OUTDIR"]
query_label = os.environ["QUERY_LABEL"]
target_label = os.environ["TARGET_LABEL"]

query_gtf = os.environ["QUERY_GTF"]
target_gtf = os.environ["TARGET_GTF"]
sam = os.environ["SAM"]
use_query_gtf = os.environ["USE_QUERY_GTF"] == "1"

min_mapq = int(os.environ["MIN_MAPQ"])
min_qcov = float(os.environ["MIN_QCOV"])
min_exon_overlap_bases = int(os.environ["MIN_EXON_OVERLAP_BASES"])
min_exon_overlap_frac = float(os.environ["MIN_EXON_OVERLAP_FRAC"])
min_best_margin = float(os.environ["MIN_BEST_MARGIN"])

cigar_re = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_attrs(attr):
    d = {}
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', attr):
        d[m.group(1)] = m.group(2)
    for part in attr.strip().split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            d.setdefault(k.strip(), v.strip())
    return d


def transcript_to_gene(gtf):
    tx2gene = {}
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = parse_attrs(fields[8])
            gene = attrs.get("gene_id") or attrs.get("gene") or attrs.get("Name")
            tx = attrs.get("transcript_id") or attrs.get("transcript") or attrs.get("ID")
            if tx and gene:
                tx2gene[tx] = gene
    return tx2gene


def infer_query_gene(transcript_id):
    m = re.match(r"(.+_g\d+)_i\d+(?:\..*)?$", transcript_id)
    if m:
        return m.group(1)
    m = re.match(r"(.+?)(?:[._-]t\d+)$", transcript_id)
    if m:
        return m.group(1)
    return transcript_id


def parse_target_annotation(gtf):
    exons = defaultdict(list)
    gene_spans = defaultdict(list)
    gene_names = {}

    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attr = fields
            attrs = parse_attrs(attr)
            gene_id = attrs.get("gene_id") or attrs.get("gene") or attrs.get("Name")
            if not gene_id:
                continue
            gene_names.setdefault(gene_id, attrs.get("gene", gene_id))
            start0 = int(start) - 1
            end0 = int(end)
            if feature == "exon":
                exons[chrom].append((start0, end0, strand, gene_id))
            elif feature == "gene":
                gene_spans[chrom].append((start0, end0, strand, gene_id))

    return make_index(exons), make_index(gene_spans), gene_names


def make_index(intervals_by_chrom):
    index = {}
    for chrom, intervals in intervals_by_chrom.items():
        intervals.sort(key=lambda x: x[0])
        starts = [x[0] for x in intervals]
        max_ends = []
        max_end = 0
        for interval in intervals:
            max_end = max(max_end, interval[1])
            max_ends.append(max_end)
        index[chrom] = (starts, intervals, max_ends)
    return index


def overlap_index(index, chrom, start, end, strand=None):
    if chrom not in index:
        return []
    starts, intervals, max_ends = index[chrom]
    hits = []
    i = bisect.bisect_left(starts, end)
    j = i - 1
    while j >= 0 and max_ends[j] > start:
        iv_start, iv_end, iv_strand, gene_id = intervals[j]
        if strand is None or iv_strand == "." or iv_strand == strand:
            ov = min(end, iv_end) - max(start, iv_start)
            if ov > 0:
                hits.append((gene_id, ov))
        j -= 1
    return hits


def cigar_blocks(pos0, cigar):
    ref_pos = pos0
    blocks = []
    q_aligned = 0
    ref_aligned = 0

    for length_s, op in cigar_re.findall(cigar):
        length = int(length_s)
        if op in ("M", "=", "X"):
            blocks.append((ref_pos, ref_pos + length))
            ref_pos += length
            q_aligned += length
            ref_aligned += length
        elif op == "I":
            q_aligned += length
        elif op == "D":
            ref_pos += length
            ref_aligned += length
        elif op == "N":
            ref_pos += length
        elif op == "S":
            pass
        elif op in ("H", "P"):
            pass

    return blocks, q_aligned, ref_aligned


def sam_records(path):
    with open(path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            seq = fields[9]
            if flag & 4 or flag & 256 or flag & 2048 or rname == "*" or cigar == "*":
                continue
            strand = "-" if flag & 16 else "+"
            qlen = len(seq)
            yield qname, rname, pos - 1, mapq, cigar, strand, qlen


query_tx2gene = transcript_to_gene(query_gtf) if use_query_gtf else {}
target_exons, target_gene_spans, target_gene_names = parse_target_annotation(target_gtf)

tx_assignments = []
unassigned = []

for qname, chrom, start0, mapq, cigar, strand, qlen in sam_records(sam):
    blocks, q_aligned, ref_aligned = cigar_blocks(start0, cigar)
    if qlen == 0:
        continue
    qcov = q_aligned / qlen
    if mapq < min_mapq or qcov < min_qcov:
        unassigned.append((qname, query_tx2gene.get(qname, infer_query_gene(qname)), chrom, start0 + 1, mapq, qcov, "low_alignment_confidence"))
        continue

    exon_overlap_by_gene = defaultdict(int)
    span_overlap_by_gene = defaultdict(int)
    for block_start, block_end in blocks:
        for gene_id, ov in overlap_index(target_exons, chrom, block_start, block_end, strand):
            exon_overlap_by_gene[gene_id] += ov
        for gene_id, ov in overlap_index(target_gene_spans, chrom, block_start, block_end, strand):
            span_overlap_by_gene[gene_id] += ov

    candidates = sorted(
        exon_overlap_by_gene.items(),
        key=lambda item: (item[1], span_overlap_by_gene.get(item[0], 0), item[0]),
        reverse=True,
    )

    if not candidates:
        span_candidates = sorted(
            span_overlap_by_gene.items(),
            key=lambda item: (item[1], item[0]),
            reverse=True,
        )
        if not span_candidates:
            unassigned.append((qname, query_tx2gene.get(qname, infer_query_gene(qname)), chrom, start0 + 1, mapq, qcov, "no_target_annotation_overlap"))
            continue

        best_gene, best_span_overlap = span_candidates[0]
        second_span_overlap = span_candidates[1][1] if len(span_candidates) > 1 else 0
        best_margin = (
            (best_span_overlap - second_span_overlap) / best_span_overlap
            if best_span_overlap
            else 0.0
        )
        query_gene = query_tx2gene.get(qname, infer_query_gene(qname))
        tx_assignments.append({
            "query_gene": query_gene,
            "query_transcript": qname,
            "target_gene": best_gene,
            "target_gene_name": target_gene_names.get(best_gene, best_gene),
            "status": "gene_span_only",
            "chrom": chrom,
            "alignment_start": start0 + 1,
            "strand": strand,
            "mapq": mapq,
            "query_coverage": qcov,
            "aligned_query_bases": q_aligned,
            "aligned_ref_bases": ref_aligned,
            "best_exon_overlap_bases": 0,
            "best_exon_overlap_frac": 0.0,
            "second_exon_overlap_bases": 0,
            "best_margin": best_margin,
            "n_overlapping_target_genes": len(span_candidates),
            "overlapping_target_genes": ",".join(f"{g}:span:{ov}" for g, ov in span_candidates[:10]),
        })
        continue

    best_gene, best_exon_overlap = candidates[0]
    second_exon_overlap = candidates[1][1] if len(candidates) > 1 else 0
    exon_overlap_frac = best_exon_overlap / q_aligned if q_aligned else 0.0
    best_margin = (
        (best_exon_overlap - second_exon_overlap) / best_exon_overlap
        if best_exon_overlap
        else 0.0
    )

    if best_exon_overlap < min_exon_overlap_bases or exon_overlap_frac < min_exon_overlap_frac:
        status = "weak_exon_overlap"
    elif second_exon_overlap and best_margin < min_best_margin:
        status = "ambiguous_exon_overlap"
    else:
        status = "assigned"

    query_gene = query_tx2gene.get(qname, infer_query_gene(qname))
    tx_assignments.append({
        "query_gene": query_gene,
        "query_transcript": qname,
        "target_gene": best_gene,
        "target_gene_name": target_gene_names.get(best_gene, best_gene),
        "status": status,
        "chrom": chrom,
        "alignment_start": start0 + 1,
        "strand": strand,
        "mapq": mapq,
        "query_coverage": qcov,
        "aligned_query_bases": q_aligned,
        "aligned_ref_bases": ref_aligned,
        "best_exon_overlap_bases": best_exon_overlap,
        "best_exon_overlap_frac": exon_overlap_frac,
        "second_exon_overlap_bases": second_exon_overlap,
        "best_margin": best_margin,
        "n_overlapping_target_genes": len(candidates),
        "overlapping_target_genes": ",".join(f"{g}:{ov}" for g, ov in candidates[:10]),
    })


gene_pair_stats = defaultdict(lambda: {
    "n_assigned_transcripts": 0,
    "n_supporting_transcripts": 0,
    "best_query_coverage": 0.0,
    "best_exon_overlap_frac": 0.0,
    "sum_exon_overlap_bases": 0,
    "sum_aligned_query_bases": 0,
    "examples": [],
    "statuses": defaultdict(int),
})

for row in tx_assignments:
    key = (row["query_gene"], row["target_gene"])
    s = gene_pair_stats[key]
    if row["status"] == "assigned":
        s["n_assigned_transcripts"] += 1
    if row["status"] in ("assigned", "ambiguous_exon_overlap", "weak_exon_overlap", "gene_span_only"):
        s["n_supporting_transcripts"] += 1
    s["best_query_coverage"] = max(s["best_query_coverage"], row["query_coverage"])
    s["best_exon_overlap_frac"] = max(s["best_exon_overlap_frac"], row["best_exon_overlap_frac"])
    s["sum_exon_overlap_bases"] += row["best_exon_overlap_bases"]
    s["sum_aligned_query_bases"] += row["aligned_query_bases"]
    s["statuses"][row["status"]] += 1
    if len(s["examples"]) < 5:
        s["examples"].append(row["query_transcript"])


best_gene_for_query = {}
for (query_gene, target_gene), s in gene_pair_stats.items():
    rank = (
        s["n_assigned_transcripts"],
        s["n_supporting_transcripts"],
        s["sum_exon_overlap_bases"],
        s["best_exon_overlap_frac"],
        s["best_query_coverage"],
    )
    if query_gene not in best_gene_for_query or rank > best_gene_for_query[query_gene][0]:
        best_gene_for_query[query_gene] = (rank, target_gene, s)


with open(f"{outdir}/query_transcript_to_target_gene_overlap.tsv", "w") as out:
    columns = [
        "query_gene",
        "query_transcript",
        "target_gene",
        "target_gene_name",
        "status",
        "chrom",
        "alignment_start",
        "strand",
        "mapq",
        "query_coverage",
        "aligned_query_bases",
        "aligned_ref_bases",
        "best_exon_overlap_bases",
        "best_exon_overlap_frac",
        "second_exon_overlap_bases",
        "best_margin",
        "n_overlapping_target_genes",
        "overlapping_target_genes",
    ]
    out.write("\t".join(columns) + "\n")
    for row in sorted(tx_assignments, key=lambda r: (r["query_gene"], r["query_transcript"], r["target_gene"])):
        out.write("\t".join(str(row[c]) if not isinstance(row[c], float) else f"{row[c]:.4f}" for c in columns) + "\n")


with open(f"{outdir}/query_to_target_gene_map.tsv", "w") as out:
    out.write(
        "query_gene\ttarget_gene\tn_assigned_transcripts\tn_supporting_transcripts\t"
        "best_query_coverage\tbest_exon_overlap_frac\tsum_exon_overlap_bases\t"
        "sum_aligned_query_bases\tstatus_counts\texample_query_transcripts\n"
    )
    for query_gene in sorted(best_gene_for_query):
        _, target_gene, s = best_gene_for_query[query_gene]
        status_counts = ",".join(f"{k}:{v}" for k, v in sorted(s["statuses"].items()))
        out.write(
            f"{query_gene}\t{target_gene}\t{s['n_assigned_transcripts']}\t"
            f"{s['n_supporting_transcripts']}\t{s['best_query_coverage']:.4f}\t"
            f"{s['best_exon_overlap_frac']:.4f}\t{s['sum_exon_overlap_bases']}\t"
            f"{s['sum_aligned_query_bases']}\t{status_counts}\t{','.join(s['examples'])}\n"
        )


with open(f"{outdir}/unassigned_query_transcripts.tsv", "w") as out:
    out.write("query_transcript\tquery_gene\tchrom\talignment_start\tmapq\tquery_coverage\treason\n")
    for row in sorted(unassigned):
        qname, query_gene, chrom, start, mapq, qcov, reason = row
        out.write(f"{qname}\t{query_gene}\t{chrom}\t{start}\t{mapq}\t{qcov:.4f}\t{reason}\n")


print(f"Wrote {outdir}/query_to_target_gene_map.tsv")
print(f"Wrote {outdir}/query_transcript_to_target_gene_overlap.tsv")
print(f"Wrote {outdir}/unassigned_query_transcripts.tsv")
PY

echo "Done."
