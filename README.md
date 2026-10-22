# KOMB: Ty1 retrotransposon insertion analysis in *S. cerevisiae*

## Purpose

This project maps and quantifies genome-wide insertion sites of the yeast
**Ty1 retrotransposon** (*Saccharomyces cerevisiae*, reference S288C) and asks
how insertion is influenced by **RNA polymerase III (Pol III)** and its
subunits. Ty1 preferentially integrates upstream of Pol III-transcribed genes
(notably tRNA genes). The analysis compares libraries from wild-type and
Pol III subunit mutants (**RPC160, RPC53, RPC34**) grown on different sugars (galactose, which
supports Ty1 induction/growth, versus dextrose, which does not) to evaluate:

- Ty1 insertion at tRNA genes: butterfly (strand-resolved, around-gene) plots,
  strand-bias box plots, per-gene and per-anticodon heatmaps, and
  mutant-vs-wild-type t-tests.
- Insertion at other Pol III-transcribed genes and at Pol II-transcribed genes
  (including nucleosome-position context).
- Genome-wide insertion profiles (heatmaps, bar plots, TPM vs coverage), the
  rDNA locus, and genes located near existing Ty1 elements.
- Results for both all alignments and uniquely mapping alignments (MAPQ >= 30).

The sole deliverable is the rendered report `Ty1_insertion_report.Rmd`.
A manuscript describing this work, including library preparation and methods,
is forthcoming; please cite it once published.

## Data

Raw sequencing data are in the NCBI SRA under BioProject
**[PRJNA1530626](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1530626)**.
Sample metadata are in `SRA_metadata.xlsx`.

Download the FASTQ files (e.g. with `sra-tools`' `prefetch`/`fasterq-dump`,
then gzip) and place them in `data/input/FASTQ/` named
`<sample>_100bp.concat_chastity_passed*.fastq.gz`; the sample name is the
part before `_100bp...` (see `data/input/hiseq_sample_id_map.tsv`).

Other required files under `data/input/` include: `genome/` (BWA index with
prefix `reference_genome`), `trna_regions.edited.bed`, `filter_regions.bed`,
`buffer_regions.tsv`, `chromosome_lengths.tsv`, `custom_feature_annotations.xlsx`,
`library_size_data.xlsx`, `hiseq_sample_id_map.tsv`, the nucleosome data
`13059_2018_1398_MOESM2_ESM.xlsx`, and `BAMs/` (extra BAMs merged with new ones).
The S288C GFF annotation is downloaded by the pipeline.

## Installation

All required software (bioinformatics tools, Python and R packages) is listed in
`environment.yml`.

```bash
conda env create -f environment.yml
conda activate komb
```

Network access is needed for the annotation download step.

## Running the pipeline

From the repository root:

```bash
./run.sh
```

`run.sh` uses `data/` as the working directory and 4 CPUs (`N_CPUS`; edit the
top of the script to change). Steps:

1. **Align** (`align_reads_to_reference_genome.sh`): trim the 5' `AACA` and the
   3' adapter with cutadapt, map with `bwa mem`, keep primary alignments with
   NM <= 1, mark/remove duplicates (Picard), and make a MAPQ >= 30
   "unique" BAM.
2. **Merge and filter BAMs** (`merge_bam_files.sh`, `filter_bam_files.sh`),
   removing regions in `filter_regions.bed`.
3. **Read count summary** (`summarise_read_counts_into_table.sh`) ->
   `data/output/summarised_read_counts.tsv`.
4. **Genome features**: download S288C annotation, convert the GFF to BED with
   buffers, define Pol II/Pol III feature sets, and find genes near Ty1.
5. **Coverage** of all features, Pol III features, Pol II positions, and the
   rDNA locus (`data/output/genome_coverage/`).
6. **Render the report** with `rmarkdown::render('./Ty1_insertion_report.Rmd')`.

Outputs go to `data/output/`; the report is written next to the Rmd and
figures to `Ty1_insertion_report_files/`.

To re-render only the report once coverage files exist:

```bash
Rscript -e "rmarkdown::render('Ty1_insertion_report.Rmd')"
```

## Repository layout

- `run.sh` - end-to-end pipeline
- `scripts/` - shell/Python helpers (`scripts/R` for R code)
- `Ty1_insertion_report.Rmd` - analysis and figures
- `environment.yml` - conda environment
- `data/input`, `data/output` - inputs and results

## Notes and known caveats

- `run.sh` creates symlinks (`data/input/alignments`, and temporary ones in
  `genome_coverage`); re-running fails on `ln -s` if they already exist.
- The `input_dir` existence check in `run.sh` occurs before the variable is
  set, so it does not actually validate anything.
- The script requires `bedtools` on the `PATH`.
