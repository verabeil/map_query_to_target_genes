# Query-to-Target Gene Mapping by Genome Overlap
Vibe coded with Codex/ChatGPT

This workflow maps genes from one annotation to another by aligning query transcripts to a target genome, then assigning each query transcript to the target gene whose annotated exons overlap the alignment the most.

The script is:

```bash
map_v2_to_ncbi_by_genome_overlap.sh
```

Despite the filename, the script is generic. It can compare V2 to NCBI, V3 to V4, or any other query/target pair with genomes and annotations.

## Requirements

Install these in the conda environment used by the script:

```bash
conda create -n v2_ncbi_map -c bioconda -c conda-forge gffread minimap2 python
```

The script currently activates this environment by default:

```bash
source /home/vbeilins/miniconda3/bin/activate
conda activate v2_ncbi_map
```

If another user has a different conda install or environment name, pass `--conda-init` and `--conda-env`, or use `--no-conda` if the tools are already on `PATH`.

## Inputs

For any comparison, you need four files:

```text
query genome FASTA, optionally .gz
query annotation GTF/GFF
target genome FASTA, optionally .gz
target annotation GTF/GFF
```

The script expects the query annotation to contain transcript IDs and gene IDs. The target annotation should contain `gene` and preferably `exon` features.

If the query is a transcriptome assembly instead of a genome annotation, use `--query-transcripts`. In that mode, the query genome and query annotation are not required.

## Current Default Run: V2 to NCBI

The defaults are set for this comparison:

```text
query genome:      Lachesis_assembly.fasta.gz
query annotation:  eupsc_models_v2.2.gtf
target genome:     GCF_053919585.1_ASM5391958v1_genomic.fna.gz
target annotation: ncbi.gtf
```

To run on the HPC:

```bash
sbatch map_v2_to_ncbi_by_genome_overlap.sh
```

Outputs will be written to:

```text
v2_to_ncbi_genome_overlap/
```

## New Comparison Example: V3 to V4

Put the script and input files in the same directory, then submit:

```bash
sbatch map_v2_to_ncbi_by_genome_overlap.sh \
  --query-label v3 \
  --query-genome V3.fasta.gz \
  --query-gtf V3.gtf \
  --target-label v4 \
  --target-genome V4.fasta.gz \
  --target-gtf V4.gtf \
  --outdir v3_to_v4_genome_overlap
```

If the conda environment has a different name:

```bash
sbatch map_v2_to_ncbi_by_genome_overlap.sh \
  --query-label v3 \
  --query-genome V3.fasta.gz \
  --query-gtf V3.gtf \
  --target-label v4 \
  --target-genome V4.fasta.gz \
  --target-gtf V4.gtf \
  --outdir v3_to_v4_genome_overlap \
  --conda-env my_mapping_env
```

## Transcriptome Example: V1 Trinity Transcriptome to V3

If V1 is only a transcriptome FASTA from Trinity, use it directly as the query:

```bash
sbatch map_v2_to_ncbi_by_genome_overlap.sh \
  --query-label v1 \
  --query-transcripts Trinity.fasta \
  --target-label v3 \
  --target-genome V3.fasta.gz \
  --target-gtf V3.gtf \
  --outdir v1_transcriptome_to_v3_genome_overlap
```

The transcriptome FASTA can also be gzipped:

```bash
--query-transcripts Trinity.fasta.gz
```

For Trinity IDs such as:

```text
TRINITY_DN123_c0_g1_i1
TRINITY_DN123_c0_g1_i2
```

the script collapses isoforms to the inferred query gene:

```text
TRINITY_DN123_c0_g1
```

If the transcriptome has non-Trinity IDs and no query annotation, each transcript ID is treated as its own query gene unless it ends with a common isoform suffix like `.t1` or `_t1`.

## Outputs

The main output is:

```text
query_to_target_gene_map.tsv
```

This contains one best target gene per query gene.

Detailed transcript-level assignments are in:

```text
query_transcript_to_target_gene_overlap.tsv
```

Transcripts that aligned poorly or did not overlap the target annotation are in:

```text
unassigned_query_transcripts.tsv
```

## Important Columns

In `query_to_target_gene_map.tsv`:

```text
query_gene
target_gene
n_assigned_transcripts
n_supporting_transcripts
best_query_coverage
best_exon_overlap_frac
sum_exon_overlap_bases
status_counts
example_query_transcripts
```

In `query_transcript_to_target_gene_overlap.tsv`, the `status` column is especially useful:

```text
assigned                 clear exon overlap to a target gene
ambiguous_exon_overlap   more than one target gene overlaps similarly
weak_exon_overlap        exon overlap exists but is below threshold
gene_span_only           overlaps a target gene span but not target exons
low_alignment_confidence MAPQ or query coverage failed
no_target_annotation_overlap good alignment but no target gene/exon overlap
```

## Thresholds

Defaults:

```bash
--min-mapq 20
--min-qcov 0.50
--min-exon-bases 50
--min-exon-frac 0.10
--min-best-margin 0.10
```

Example with stricter mapping:

```bash
sbatch map_v2_to_ncbi_by_genome_overlap.sh \
  --query-label v3 \
  --query-genome V3.fasta.gz \
  --query-gtf V3.gtf \
  --target-label v4 \
  --target-genome V4.fasta.gz \
  --target-gtf V4.gtf \
  --outdir v3_to_v4_strict_overlap \
  --min-mapq 30 \
  --min-qcov 0.80 \
  --min-exon-frac 0.30
```

## Checking a Job

Submit:

```bash
sbatch map_v2_to_ncbi_by_genome_overlap.sh
```

Check status:

```bash
squeue -u vbeilins
```

Check the SLURM log files after completion:

```bash
ls -lh map_v2_ncbi_loci.out map_v2_ncbi_loci.err
```

If the job fails early, first check that all input filenames match exactly and that `gffread`, `minimap2`, and `python3` are available in the activated environment.
